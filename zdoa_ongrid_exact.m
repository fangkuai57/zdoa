function [spectrum, doa_est, debug] = zdoa_ongrid_exact(z, W, scan_deg, num_sources, cfg)
% Exact on-grid DOA recovery from averaged beam powers z in R^{P x 1}.
% This version is intended for simple on-grid scenes with small K.

[P, M] = size(W);
N = numel(scan_deg);
L_snap = get_cfg(cfg, 'L', 1);
max_support_evals = get_cfg(cfg, 'max_support_evals', 50000);
candidate_pool_size = get_cfg(cfg, 'candidate_pool_size', max(16, 4 * num_sources));

D = zeros(P, N);
for n = 1:N
    a = steering_vector_local(scan_deg(n), M, cfg.d_over_lambda);
    D(:, n) = abs(W * a).^2;
end
c = sum(abs(W).^2, 2);

num_combos = nchoosek_safe(N, num_sources);
if num_combos <= max_support_evals
    candidate_idx = 1:N;
else
    matched = max(D' * z, 0);
    [~, order] = sort(matched, 'descend');
    candidate_idx = sort(order(1:min(candidate_pool_size, numel(order))));
end

combos = nchoosek(candidate_idx, num_sources);
best_score = inf;
best_idx = combos(1, :);
best_p = [];
best_sigma2 = 0;
best_u = [];

for combo_idx = 1:size(combos, 1)
    support_idx = combos(combo_idx, :);
    [p_support, sigma2_est, u_model, score] = fit_support_gamma(z, D(:, support_idx), c, L_snap, cfg);
    if score < best_score
        best_score = score;
        best_idx = support_idx;
        best_p = p_support;
        best_sigma2 = sigma2_est;
        best_u = u_model;
    end
end

spectrum = zeros(N, 1);
spectrum(best_idx) = best_p;
if max(spectrum) > 0
    spectrum = spectrum / max(spectrum);
end

doa_est = sort(scan_deg(best_idx));
debug.D = D;
debug.c = c;
debug.best_support_idx = best_idx;
debug.best_support_deg = scan_deg(best_idx);
debug.best_p = best_p;
debug.best_sigma2 = best_sigma2;
debug.best_u = best_u;
debug.best_score = best_score;
debug.num_evaluated_supports = size(combos, 1);
debug.frontend_output_size = [P, 1];
end

function [p_support, sigma2_est, u_model, score] = fit_support_gamma(z, Ds, c, L_snap, cfg)
max_alt_iters = get_cfg(cfg, 'max_alt_iters', 12);
fit_floor = get_cfg(cfg, 'fit_floor', 1e-8);

x = lsqnonneg([Ds, c], z);
p_support = x(1:end-1);
sigma2_est = x(end);

for iter = 1:max_alt_iters
    u_model = max(Ds * p_support + sigma2_est * c, fit_floor);
    weights = 1 ./ max(u_model .^ 2, fit_floor);
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

u_model = max(Ds * p_support + sigma2_est * c, fit_floor);
score = sum(L_snap * log(u_model) + L_snap * (z ./ u_model));
end

function a = steering_vector_local(theta_deg, M, d_over_lambda)
n = (0:M-1).';
a = exp(1i * 2 * pi * d_over_lambda * n * sind(theta_deg));
a = a / sqrt(M);
end

function value = get_cfg(cfg, field_name, default_value)
if isfield(cfg, field_name)
    value = cfg.(field_name);
else
    value = default_value;
end
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
