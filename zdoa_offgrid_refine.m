function [spectrum, doa_est, debug] = zdoa_offgrid_refine(z, W, scan_deg, num_sources, cfg)
% Power-domain off-grid DOA recovery.
% Step 1: on-grid exact support initialization.
% Step 2: coordinate-wise continuous angle refinement.

[P, M] = size(W);
cfg = fill_defaults(cfg, scan_deg);

[~, doa_init, init_debug] = zdoa_ongrid_exact(z, W, scan_deg, num_sources, cfg);
theta = sort(doa_init(:));

[theta, p_support, sigma2_est, best_score] = refine_support(theta, z, W, cfg, M);

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
debug.best_sigma2 = sigma2_est;
debug.best_score = best_score;
debug.frontend_output_size = [P, 1];
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
end

function [theta, p_support, sigma2_est, best_score] = refine_support(theta_init, z, W, cfg, M)
theta = sort(theta_init(:));
[p_support, sigma2_est, best_score] = fit_support_theta(z, W, theta, cfg, M);

for pass_idx = 1:cfg.max_refine_passes
    improved = false;
    for radius = cfg.refine_radii_deg
        for k = 1:numel(theta)
            candidate_deg = linspace(theta(k) - radius, theta(k) + radius, cfg.refine_samples);
            candidate_deg = unique([candidate_deg, theta(k)]);
            candidate_deg = candidate_deg(candidate_deg >= cfg.angle_min_deg & candidate_deg <= cfg.angle_max_deg);

            best_local_theta = theta;
            best_local_p = p_support;
            best_local_sigma2 = sigma2_est;
            best_local_score = best_score;

            for cand_idx = 1:numel(candidate_deg)
                theta_trial = theta;
                theta_trial(k) = candidate_deg(cand_idx);
                if violates_separation(theta_trial, k, cfg.refine_guard_deg)
                    continue;
                end
                theta_trial = sort(theta_trial);
                [p_trial, sigma2_trial, score_trial] = fit_support_theta(z, W, theta_trial, cfg, M);
                if score_trial + 1e-10 < best_local_score
                    best_local_theta = theta_trial;
                    best_local_p = p_trial;
                    best_local_sigma2 = sigma2_trial;
                    best_local_score = score_trial;
                end
            end

            if best_local_score + 1e-10 < best_score
                theta = best_local_theta;
                p_support = best_local_p;
                sigma2_est = best_local_sigma2;
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

function [p_support, sigma2_est, score] = fit_support_theta(z, W, theta, cfg, M)
Ds = zeros(size(W, 1), numel(theta));
for k = 1:numel(theta)
    a = steering_vector_local(theta(k), M, cfg.d_over_lambda);
    Ds(:, k) = abs(W * a) .^ 2;
end
c = sum(abs(W) .^ 2, 2);

x = lsqnonneg([Ds, c], z);
p_support = x(1:end-1);
sigma2_est = x(end);

for iter = 1:cfg.max_alt_iters
    u_model = max(Ds * p_support + sigma2_est * c, cfg.fit_floor);
    weights = 1 ./ max(u_model .^ 2, cfg.fit_floor);
    A = [bsxfun(@times, Ds, sqrt(weights)), sqrt(weights) .* c];
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

u_model = max(Ds * p_support + sigma2_est * c, cfg.fit_floor);
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
