function [spectrum, doa_est, debug] = teacher_codebook_decode_wideband(z_bands, W_bands, coarse_scan_deg, num_sources, cfg)
% Wideband teacher-style decoder with shared dense support.

[P, M, num_bands] = size(W_bands);
cfg = fill_defaults(cfg, coarse_scan_deg, P);

fine_grid_deg = cfg.fine_grid_deg(:).';
D_bands = zeros(P, numel(fine_grid_deg), num_bands);
c_bands = zeros(P, num_bands);
for band_idx = 1:num_bands
    d_eff = cfg.d_over_lambda * cfg.freq_ratios(band_idx);
    D_bands(:, :, band_idx) = build_power_dictionary(W_bands(:, :, band_idx), fine_grid_deg, d_eff, M);
    c_bands(:, band_idx) = sum(abs(W_bands(:, :, band_idx)) .^ 2, 2);
end

[candidate_idx, single_scores] = build_candidate_pool(z_bands, D_bands, c_bands, num_sources, cfg);
[best_idx, best_p, best_sigma, best_score] = exact_support_decode(z_bands, D_bands, c_bands, candidate_idx, num_sources, cfg);
[best_theta, best_p, best_sigma, best_score] = local_grid_refine(z_bands, W_bands, fine_grid_deg(best_idx), best_p, best_sigma, best_score, cfg, M);

doa_est = sort(best_theta(:)).';
spectrum = zeros(numel(coarse_scan_deg), 1);
for k = 1:numel(doa_est)
    [~, idx] = min(abs(coarse_scan_deg - doa_est(k)));
    spectrum(idx) = spectrum(idx) + best_p(k);
end
peak = max(spectrum);
if peak > 0
    spectrum = spectrum / peak;
end

debug.fine_grid_deg = fine_grid_deg;
debug.single_scores = single_scores;
debug.candidate_idx = candidate_idx;
debug.best_theta = doa_est;
debug.best_p = best_p(:).';
debug.best_sigma = best_sigma(:).';
debug.best_score = best_score;
debug.frontend_output_size = [P, 1];
debug.num_bands = num_bands;
end

function cfg = fill_defaults(cfg, coarse_scan_deg, P)
cfg = set_default(cfg, 'd_over_lambda', 0.5);
cfg = set_default(cfg, 'L', 40);
cfg = set_default(cfg, 'fit_floor', 1e-8);
cfg = set_default(cfg, 'max_alt_iters', 12);
cfg = set_default(cfg, 'fine_step_deg', 0.05);
cfg = set_default(cfg, 'pool_size', max(12, 4 * P));
cfg = set_default(cfg, 'pool_guard_deg', 0.35);
cfg = set_default(cfg, 'refine_radii_deg', [0.30, 0.15, 0.08, 0.04]);
cfg = set_default(cfg, 'refine_samples', 9);
cfg = set_default(cfg, 'refine_guard_deg', 0.25);
cfg = set_default(cfg, 'angle_min_deg', min(coarse_scan_deg));
cfg = set_default(cfg, 'angle_max_deg', max(coarse_scan_deg));

if ~isfield(cfg, 'fine_grid_deg') || isempty(cfg.fine_grid_deg)
    cfg.fine_grid_deg = cfg.angle_min_deg:cfg.fine_step_deg:cfg.angle_max_deg;
end
end

function D = build_power_dictionary(W, grid_deg, d_over_lambda, M)
D = zeros(size(W, 1), numel(grid_deg));
for g = 1:numel(grid_deg)
    a = steering_vector_local(grid_deg(g), M, d_over_lambda);
    D(:, g) = abs(W * a) .^ 2;
end
end

function [candidate_idx, single_scores] = build_candidate_pool(z_bands, D_bands, c_bands, num_sources, cfg)
num_grid = size(D_bands, 2);
num_bands = size(D_bands, 3);
single_scores = inf(num_grid, 1);
for g = 1:num_grid
    [~, ~, score] = fit_support(z_bands, D_bands(:, g, :), c_bands, cfg);
    single_scores(g) = score;
end

[~, order] = sort(single_scores, 'ascend');
candidate_idx = zeros(1, 0);
for idx = order(:).'
    theta = cfg.fine_grid_deg(idx);
    if isempty(candidate_idx) || all(abs(theta - cfg.fine_grid_deg(candidate_idx)) >= cfg.pool_guard_deg)
        candidate_idx(end + 1) = idx; %#ok<AGROW>
    end
    if numel(candidate_idx) >= max(cfg.pool_size, 4 * num_sources + num_bands)
        break;
    end
end

if numel(candidate_idx) < num_sources
    candidate_idx = sort(order(1:num_sources));
else
    candidate_idx = sort(candidate_idx);
end
end

function [best_idx, best_p, best_sigma, best_score] = exact_support_decode(z_bands, D_bands, c_bands, candidate_idx, num_sources, cfg)
combos = nchoosek(candidate_idx, num_sources);
best_score = inf;
best_idx = combos(1, :);
best_p = zeros(num_sources, 1);
best_sigma = zeros(size(z_bands, 2), 1);

for row_idx = 1:size(combos, 1)
    support_idx = combos(row_idx, :);
    [p_support, sigma_vec, score] = fit_support(z_bands, D_bands(:, support_idx, :), c_bands, cfg);
    if score < best_score
        best_score = score;
        best_idx = support_idx;
        best_p = p_support;
        best_sigma = sigma_vec;
    end
end
end

function [theta, p_support, sigma_vec, best_score] = local_grid_refine(z_bands, W_bands, theta_init, p_support, sigma_vec, best_score, cfg, M)
theta = sort(theta_init(:));
for radius = cfg.refine_radii_deg
    improved = false;
    for k = 1:numel(theta)
        candidate_deg = linspace(theta(k) - radius, theta(k) + radius, cfg.refine_samples);
        candidate_deg = unique([candidate_deg, theta(k)]);
        candidate_deg = candidate_deg(candidate_deg >= cfg.angle_min_deg & candidate_deg <= cfg.angle_max_deg);

        for cand = candidate_deg
            theta_trial = theta;
            theta_trial(k) = cand;
            if violates_separation(theta_trial, k, cfg.refine_guard_deg)
                continue;
            end
            theta_trial = sort(theta_trial);
            D_trial = zeros(size(W_bands, 1), numel(theta_trial), size(W_bands, 3));
            c_trial = zeros(size(W_bands, 1), size(W_bands, 3));
            for band_idx = 1:size(W_bands, 3)
                d_eff = cfg.d_over_lambda * cfg.freq_ratios(band_idx);
                for j = 1:numel(theta_trial)
                    a = steering_vector_local(theta_trial(j), M, d_eff);
                    D_trial(:, j, band_idx) = abs(W_bands(:, :, band_idx) * a) .^ 2;
                end
                c_trial(:, band_idx) = sum(abs(W_bands(:, :, band_idx)) .^ 2, 2);
            end
            [p_trial, sigma_trial, score_trial] = fit_support(z_bands, D_trial, c_trial, cfg);
            if score_trial + 1e-10 < best_score
                theta = theta_trial;
                p_support = p_trial;
                sigma_vec = sigma_trial;
                best_score = score_trial;
                improved = true;
            end
        end
    end
    if ~improved
        continue;
    end
end
end

function tf = violates_separation(theta, idx, guard_deg)
others = theta;
others(idx) = [];
tf = ~isempty(others) && any(abs(theta(idx) - others) < guard_deg);
end

function [p_support, sigma_vec, score] = fit_support(z_bands, D_support, c_bands, cfg)
[A, y] = build_system(z_bands, D_support, c_bands);
x = lsqnonneg(A, y);
K = size(D_support, 2);
p_support = x(1:K);
sigma_vec = x(K + 1:end);

for iter = 1:cfg.max_alt_iters
    u_model = max(A * [p_support; sigma_vec], cfg.fit_floor);
    weights = 1 ./ max(u_model .^ 2, cfg.fit_floor);
    Aw = bsxfun(@times, A, sqrt(weights));
    bw = sqrt(weights) .* y;
    x_next = lsqnonneg(Aw, bw);
    if norm(x_next - [p_support; sigma_vec]) <= 1e-8 * max(1, norm([p_support; sigma_vec]))
        p_support = x_next(1:K);
        sigma_vec = x_next(K + 1:end);
        break;
    end
    p_support = x_next(1:K);
    sigma_vec = x_next(K + 1:end);
end

u_model = max(A * [p_support; sigma_vec], cfg.fit_floor);
score = sum(cfg.L * log(u_model) + cfg.L * (y ./ u_model));
end

function [A, y] = build_system(z_bands, D_support, c_bands)
[P, K, num_bands] = size(D_support);
A_main = zeros(P * num_bands, K);
C_block = zeros(P * num_bands, num_bands);
y = zeros(P * num_bands, 1);
for band_idx = 1:num_bands
    rows = (band_idx - 1) * P + (1:P);
    A_main(rows, :) = D_support(:, :, band_idx);
    C_block(rows, band_idx) = c_bands(:, band_idx);
    y(rows) = z_bands(:, band_idx);
end
A = [A_main, C_block];
end

function a = steering_vector_local(theta_deg, M, d_over_lambda)
n = (0:M - 1).';
a = exp(1i * 2 * pi * d_over_lambda * sind(theta_deg) * n) / sqrt(M);
end

function cfg = set_default(cfg, field_name, default_value)
if ~isfield(cfg, field_name)
    cfg.(field_name) = default_value;
end
end
