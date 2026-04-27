function [spectrum, doa_est, debug] = teacher_codebook_decode(z, W, coarse_scan_deg, num_sources, cfg)
% Teacher-style decoder for localized/codebook beamformers.
% It uses a dense off-grid codebook and exact support search on a pruned pool.

[P, M] = size(W);
cfg = fill_defaults(cfg, coarse_scan_deg, P);

fine_grid_deg = cfg.fine_grid_deg(:).';
D_fine = build_power_dictionary(W, fine_grid_deg, cfg.d_over_lambda, M);
c = sum(abs(W) .^ 2, 2);

[candidate_idx, single_scores] = build_candidate_pool(z, D_fine, c, num_sources, cfg);
[best_idx, best_p, best_sigma2, best_score] = exact_support_decode(z, D_fine, c, candidate_idx, num_sources, cfg);
[best_theta, best_p, best_sigma2, best_score] = local_grid_refine(z, W, fine_grid_deg(best_idx), best_p, best_sigma2, best_score, cfg, M);

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
debug.best_sigma2 = best_sigma2;
debug.best_score = best_score;
debug.frontend_output_size = [P, 1];
end

function cfg = fill_defaults(cfg, coarse_scan_deg, P)
cfg = set_default(cfg, 'd_over_lambda', 0.5);
cfg = set_default(cfg, 'L', 40);
cfg = set_default(cfg, 'fit_floor', 1e-8);
cfg = set_default(cfg, 'max_alt_iters', 12);
cfg = set_default(cfg, 'fine_step_deg', 0.05);
cfg = set_default(cfg, 'pool_size', max(10, 4 * P));
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

function [candidate_idx, single_scores] = build_candidate_pool(z, D_fine, c, num_sources, cfg)
num_grid = size(D_fine, 2);
single_scores = inf(num_grid, 1);
for g = 1:num_grid
    [~, ~, score] = fit_support(z, D_fine(:, g), c, cfg);
    single_scores(g) = score;
end

[~, order] = sort(single_scores, 'ascend');
candidate_idx = zeros(1, 0);
for idx = order(:).'
    theta = cfg.fine_grid_deg(idx);
    if isempty(candidate_idx) || all(abs(theta - cfg.fine_grid_deg(candidate_idx)) >= cfg.pool_guard_deg)
        candidate_idx(end + 1) = idx; %#ok<AGROW>
    end
    if numel(candidate_idx) >= max(cfg.pool_size, 3 * num_sources)
        break;
    end
end

if numel(candidate_idx) < num_sources
    candidate_idx = sort(order(1:num_sources));
else
    candidate_idx = sort(candidate_idx);
end
end

function [best_idx, best_p, best_sigma2, best_score] = exact_support_decode(z, D_fine, c, candidate_idx, num_sources, cfg)
combos = nchoosek(candidate_idx, num_sources);
best_score = inf;
best_idx = combos(1, :);
best_p = zeros(num_sources, 1);
best_sigma2 = 0;

for row_idx = 1:size(combos, 1)
    support_idx = combos(row_idx, :);
    [p_support, sigma2_est, score] = fit_support(z, D_fine(:, support_idx), c, cfg);
    if score < best_score
        best_score = score;
        best_idx = support_idx;
        best_p = p_support;
        best_sigma2 = sigma2_est;
    end
end
end

function [theta, p_support, sigma2_est, best_score] = local_grid_refine(z, W, theta_init, p_support, sigma2_est, best_score, cfg, M)
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
            D_trial = build_power_dictionary(W, theta_trial, cfg.d_over_lambda, M);
            [p_trial, sigma_trial, score_trial] = fit_support(z, D_trial, sum(abs(W) .^ 2, 2), cfg);
            if score_trial + 1e-10 < best_score
                theta = theta_trial;
                p_support = p_trial;
                sigma2_est = sigma_trial;
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

function [p_support, sigma2_est, score] = fit_support(z, D_support, c, cfg)
if isvector(D_support)
    D_support = D_support(:);
end
x = lsqnonneg([D_support, c], z);
p_support = x(1:end-1);
sigma2_est = x(end);

for iter = 1:cfg.max_alt_iters
    u_model = max(D_support * p_support + sigma2_est * c, cfg.fit_floor);
    weights = 1 ./ max(u_model .^ 2, cfg.fit_floor);
    A = [bsxfun(@times, D_support, sqrt(weights)), sqrt(weights) .* c];
    b = sqrt(weights) .* z;
    x_next = lsqnonneg(A, b);
    if norm(x_next - [p_support; sigma2_est]) <= 1e-8 * max(1, norm([p_support; sigma2_est]))
        p_support = x_next(1:end-1);
        sigma2_est = x_next(end);
        break;
    end
    p_support = x_next(1:end-1);
    sigma2_est = x_next(end);
end

u_model = max(D_support * p_support + sigma2_est * c, cfg.fit_floor);
score = sum(cfg.L * log(u_model) + cfg.L * (z ./ u_model));
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
