function [spectrum, doa_est, debug] = zdoa_wideband_offgrid_refine(z_bands, D_bands, c_bands, W_bands, scan_deg, num_sources, cfg)
% Wideband power-domain off-grid DOA recovery with shared continuous support.

[P, ~, num_bands] = size(D_bands);
cfg = fill_defaults(cfg, scan_deg);

[~, doa_init, init_debug] = wideband_init_support(z_bands, D_bands, c_bands, scan_deg, num_sources, cfg);
theta = sort(doa_init(:));

[theta, p_support, sigma_vec, best_score] = refine_support(theta, z_bands, W_bands, cfg);

spectrum = zeros(numel(scan_deg), 1);
for k = 1:numel(theta)
    [~, idx] = min(abs(scan_deg - theta(k)));
    spectrum(idx) = spectrum(idx) + p_support(k);
end
peak = max(spectrum);
if peak > 0
    spectrum = spectrum / peak;
end

doa_est = sort(theta(:)).';
debug.init_doa = doa_init;
debug.best_support_idx = init_debug.best_support_idx;
debug.best_support_deg = theta(:).';
debug.best_p = p_support(:).';
debug.best_sigma = sigma_vec(:).';
debug.best_score = best_score;
debug.frontend_output_size = [P, 1];
debug.num_bands = num_bands;
end

function cfg = fill_defaults(cfg, scan_deg)
cfg = set_default(cfg, 'L', 1);
cfg = set_default(cfg, 'd_over_lambda', 0.5);
cfg = set_default(cfg, 'fit_floor', 1e-8);
cfg = set_default(cfg, 'max_alt_iters', 12);
cfg = set_default(cfg, 'max_refine_passes', 2);
cfg = set_default(cfg, 'refine_radii_deg', [1.0, 0.5, 0.25, 0.12, 0.06]);
cfg = set_default(cfg, 'refine_samples', 11);
cfg = set_default(cfg, 'refine_guard_deg', 0.6);
cfg = set_default(cfg, 'angle_min_deg', min(scan_deg));
cfg = set_default(cfg, 'angle_max_deg', max(scan_deg));
cfg = set_default(cfg, 'max_support_evals', 50000);
cfg = set_default(cfg, 'candidate_pool_size', max(20, 4 *  numel(cfg.refine_radii_deg)));
end

function [spectrum, doa_init, debug] = wideband_init_support(z_bands, D_bands, c_bands, scan_deg, num_sources, cfg)
num_grid = numel(scan_deg);
num_combos = nchoosek_safe(num_grid, num_sources);
if num_combos <= cfg.max_support_evals
    candidate_idx = 1:num_grid;
else
    matched = zeros(num_grid, 1);
    for band_idx = 1:size(z_bands, 2)
        matched = matched + max(D_bands(:, :, band_idx)' * z_bands(:, band_idx), 0);
    end
    [~, order] = sort(matched, 'descend');
    pool_size = min(max(20, 4 * num_sources), num_grid);
    candidate_idx = sort(order(1:pool_size));
end

combos = nchoosek(candidate_idx, num_sources);
best_score = inf;
best_idx = combos(1, :);
best_p = zeros(num_sources, 1);
best_sigma = zeros(size(z_bands, 2), 1);

for combo_idx = 1:size(combos, 1)
    support_idx = combos(combo_idx, :);
    [p_support, sigma_vec, score] = fit_support_wideband(z_bands, D_bands(:, support_idx, :), c_bands, cfg);
    if score < best_score
        best_score = score;
        best_idx = support_idx;
        best_p = p_support;
        best_sigma = sigma_vec;
    end
end

spectrum = zeros(num_grid, 1);
spectrum(best_idx) = best_p;
peak = max(spectrum);
if peak > 0
    spectrum = spectrum / peak;
end

doa_init = sort(scan_deg(best_idx));
debug.best_support_idx = best_idx;
debug.best_sigma = best_sigma;
debug.best_score = best_score;
end

function [theta, p_support, sigma_vec, best_score] = refine_support(theta_init, z_bands, W_bands, cfg)
theta = sort(theta_init(:));
[p_support, sigma_vec, best_score] = fit_theta_wideband(theta, z_bands, W_bands, cfg);

for pass_idx = 1:cfg.max_refine_passes
    improved = false;
    for radius = cfg.refine_radii_deg
        for k = 1:numel(theta)
            candidate_deg = linspace(theta(k) - radius, theta(k) + radius, cfg.refine_samples);
            candidate_deg = unique([candidate_deg, theta(k)]);
            candidate_deg = candidate_deg(candidate_deg >= cfg.angle_min_deg & candidate_deg <= cfg.angle_max_deg);

            best_local_theta = theta;
            best_local_p = p_support;
            best_local_sigma = sigma_vec;
            best_local_score = best_score;

            for cand_idx = 1:numel(candidate_deg)
                theta_trial = theta;
                theta_trial(k) = candidate_deg(cand_idx);
                if violates_separation(theta_trial, k, cfg.refine_guard_deg)
                    continue;
                end
                theta_trial = sort(theta_trial);
                [p_trial, sigma_trial, score_trial] = fit_theta_wideband(theta_trial, z_bands, W_bands, cfg);
                if score_trial + 1e-10 < best_local_score
                    best_local_theta = theta_trial;
                    best_local_p = p_trial;
                    best_local_sigma = sigma_trial;
                    best_local_score = score_trial;
                end
            end

            if best_local_score + 1e-10 < best_score
                theta = best_local_theta;
                p_support = best_local_p;
                sigma_vec = best_local_sigma;
                best_score = best_local_score;
                improved = true;
            end
        end
    end

    if ~improved
        break;
    end
end
end

function tf = violates_separation(theta, idx, guard_deg)
others = theta;
others(idx) = [];
if isempty(others)
    tf = false;
    return;
end
tf = any(abs(theta(idx) - others) < guard_deg);
end

function [p_support, sigma_vec, score] = fit_theta_wideband(theta, z_bands, W_bands, cfg)
num_bands = size(W_bands, 3);
P = size(W_bands, 1);
K = numel(theta);
Ds_bands = zeros(P, K, num_bands);
c_bands = zeros(P, num_bands);

for band_idx = 1:num_bands
    d_eff = cfg.d_over_lambda * cfg.freq_ratios(band_idx);
    W_band = W_bands(:, :, band_idx);
    for k = 1:K
        a = steering_vector_local(theta(k), size(W_band, 2), d_eff);
        Ds_bands(:, k, band_idx) = abs(W_band * a) .^ 2;
    end
    c_bands(:, band_idx) = sum(abs(W_band) .^ 2, 2);
end

[p_support, sigma_vec, score] = fit_support_wideband(z_bands, Ds_bands, c_bands, cfg);
end

function [p_support, sigma_vec, score] = fit_support_wideband(z_bands, Ds_bands, c_bands, cfg)
[A, y] = build_support_system(z_bands, Ds_bands, c_bands);
x = lsqnonneg(A, y);
num_support = size(Ds_bands, 2);
p_support = x(1:num_support);
sigma_vec = x(num_support + 1:end);

for iter = 1:cfg.max_alt_iters
    u_model = max(A * [p_support; sigma_vec], cfg.fit_floor);
    weights = 1 ./ max(u_model .^ 2, cfg.fit_floor);
    Aw = bsxfun(@times, A, sqrt(weights));
    bw = sqrt(weights) .* y;
    x_next = lsqnonneg(Aw, bw);
    if norm(x_next - [p_support; sigma_vec]) <= 1e-8 * max(1, norm([p_support; sigma_vec]))
        p_support = x_next(1:num_support);
        sigma_vec = x_next(num_support + 1:end);
        break;
    end
    p_support = x_next(1:num_support);
    sigma_vec = x_next(num_support + 1:end);
end

u_model = max(A * [p_support; sigma_vec], cfg.fit_floor);
score = sum(cfg.L * log(u_model) + cfg.L * (y ./ u_model));
end

function [A, y] = build_support_system(z_bands, Ds_bands, c_bands)
[P, num_support, num_bands] = size(Ds_bands);
A_main = zeros(P * num_bands, num_support);
C_block = zeros(P * num_bands, num_bands);
y = zeros(P * num_bands, 1);
for band_idx = 1:num_bands
    rows = (band_idx - 1) * P + (1:P);
    A_main(rows, :) = Ds_bands(:, :, band_idx);
    C_block(rows, band_idx) = c_bands(:, band_idx);
    y(rows) = z_bands(:, band_idx);
end
A = [A_main, C_block];
end

function a = steering_vector_local(theta_deg, M, d_over_lambda)
n = (0:M - 1).';
a = exp(1i * 2 * pi * d_over_lambda * sind(theta_deg) * n) / sqrt(M);
end

function value = nchoosek_safe(n, k)
if k > n
    value = 0;
    return;
end
value = 1;
for idx = 1:k
    value = value * (n - k + idx) / idx;
    if value > 1e9
        break;
    end
end
end

function cfg = set_default(cfg, field_name, default_value)
if ~isfield(cfg, field_name)
    cfg.(field_name) = default_value;
end
end
