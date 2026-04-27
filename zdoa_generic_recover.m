function [spectrum, doa_est, debug] = zdoa_generic_recover(z, W, scan_deg, num_sources, cfg)
% Z-domain DOA recovery from averaged beam powers z in R^{P x 1}.
% The solver combines weighted sparse recovery with local support refinement.

[B, M] = size(W);
N = numel(scan_deg);

D = zeros(B, N);
for n = 1:N
    a = steering_vector_local(scan_deg(n), M, cfg.d_over_lambda);
    D(:, n) = abs(W * a).^2;
end
c = sum(abs(W).^2, 2);

[spectrum_sparse, sparse_debug] = solve_weighted_sparse(z, D, c, cfg);
peak_idx = pick_candidate_peaks(spectrum_sparse, max(get_cfg(cfg, 'candidate_pool_size', 8), num_sources));
refine_idx = refine_support(peak_idx, D, c, z, num_sources, scan_deg, cfg);
[spectrum, q_support, sigma2_est, fit_error] = fit_support(refine_idx, D, c, z, cfg);

doa_est = sort(scan_deg(refine_idx));
debug.D = D;
debug.c = c;
debug.sparse_spectrum = spectrum_sparse;
debug.sparse = sparse_debug;
debug.support_idx = refine_idx;
debug.q_support = q_support;
debug.sigma2_est = sigma2_est;
debug.fit_error = fit_error;
end

function [spectrum, debug] = solve_weighted_sparse(z, D, c, cfg)
lambda_sparse = get_cfg(cfg, 'lambda_sparse', 0.04);
lambda_smooth = get_cfg(cfg, 'lambda_smooth', 0.01);
reweight_loops = get_cfg(cfg, 'reweight_loops', 4);
reweight_epsilon = get_cfg(cfg, 'reweight_epsilon', 1e-3);
max_iters = get_cfg(cfg, 'max_iters', 350);
step_shrink = get_cfg(cfg, 'step_shrink', 0.98);
power_weight_floor = get_cfg(cfg, 'power_weight_floor', 1e-6);

N = size(D, 2);
G = diff(eye(N), 1, 1);
Lh = G' * G;

q = max(D' * z, 0);
if max(q) > 0
    q = q / max(q);
end

sparsity_weights = ones(N, 1);
measurement_weights = ones(size(z));

for outer = 1:reweight_loops
    sqrt_w = sqrt(measurement_weights);
    Dw = bsxfun(@times, D, sqrt_w);
    cw = sqrt_w .* c;
    zw = sqrt_w .* z;

    Lipschitz = norm(Dw' * Dw + lambda_smooth * Lh, 2);
    if Lipschitz <= 0
        Lipschitz = 1;
    end
    step = step_shrink / Lipschitz;

    sigma2 = max(0, (cw' * (zw - Dw * q)) / max(cw' * cw, eps));
    x = q;
    y = q;
    t = 1;

    for iter = 1:max_iters
        residual = Dw * y + sigma2 * cw - zw;
        grad = Dw' * residual + lambda_smooth * (Lh * y);
        q_next = max(0, y - step * grad - step * lambda_sparse * sparsity_weights);
        sigma2 = max(0, (cw' * (zw - Dw * q_next)) / max(cw' * cw, eps));

        t_next = 0.5 * (1 + sqrt(1 + 4 * t^2));
        y = q_next + ((t - 1) / t_next) * (q_next - x);

        if norm(q_next - x) <= 1e-6 * max(1, norm(x))
            x = q_next;
            break;
        end

        x = q_next;
        t = t_next;
    end

    q = x;
    sigma2 = max(0, (c' * (z - D * q)) / max(c' * c, eps));
    z_model = D * q + sigma2 * c;
    measurement_weights = 1 ./ max(z_model, power_weight_floor);
    measurement_weights = measurement_weights / mean(measurement_weights);
    sparsity_weights = 1 ./ (q + reweight_epsilon);
    sparsity_weights = sparsity_weights / mean(sparsity_weights);
end

spectrum = q;
if max(spectrum) > 0
    spectrum = spectrum / max(spectrum);
end

debug.q = q;
debug.sigma2_est = sigma2;
debug.measurement_weights = measurement_weights;
debug.sparsity_weights = sparsity_weights;
end

function support_idx = refine_support(candidate_idx, D, c, z, num_sources, scan_deg, cfg)
candidate_idx = unique(candidate_idx(candidate_idx > 0));
if isempty(candidate_idx)
    [~, best_idx] = max(D' * z);
    candidate_idx = best_idx;
end

scan_step = max(abs(median(diff(scan_deg))), 1e-6);
refine_radius = max(1, round(get_cfg(cfg, 'refine_span_deg', 1.0) / scan_step));
expanded_idx = candidate_idx(:).';
for idx = 1:numel(candidate_idx)
    local_idx = max(1, candidate_idx(idx) - refine_radius):min(numel(scan_deg), candidate_idx(idx) + refine_radius);
    expanded_idx = [expanded_idx, local_idx]; %#ok<AGROW>
end
expanded_idx = unique(expanded_idx);
max_refine_candidates = get_cfg(cfg, 'max_refine_candidates', 12);
if numel(expanded_idx) > max_refine_candidates
    [~, order] = sort(max(D(:, expanded_idx)' * z, 0), 'descend');
    expanded_idx = expanded_idx(order(1:max_refine_candidates));
end

if numel(expanded_idx) <= num_sources
    support_idx = sort(expanded_idx(1:min(end, num_sources)));
    return;
end

combos = nchoosek(expanded_idx, num_sources);
best_score = inf;
support_idx = combos(1, :);
for combo_idx = 1:size(combos, 1)
    [~, ~, ~, score] = fit_support(combos(combo_idx, :), D, c, z, cfg);
    if score < best_score
        best_score = score;
        support_idx = combos(combo_idx, :);
    end
end
support_idx = sort(support_idx);
end

function [spectrum, q_support, sigma2_est, score] = fit_support(support_idx, D, c, z, cfg)
N = size(D, 2);
weights = 1 ./ max(z, get_cfg(cfg, 'power_weight_floor', 1e-6));
weights = weights / mean(weights);
sqrt_w = sqrt(weights(:));
A = [bsxfun(@times, D(:, support_idx), sqrt_w), sqrt_w .* c(:)];
zz = sqrt_w .* z(:);
x = lsqnonneg(A, zz);
q_support = x(1:numel(support_idx));
sigma2_est = x(end);
residual = zz - A * x;
score = residual' * residual;
spectrum = zeros(N, 1);
spectrum(support_idx) = q_support;
if max(spectrum) > 0
    spectrum = spectrum / max(spectrum);
end
end

function idx = pick_candidate_peaks(spectrum, max_count)
s = spectrum(:);
N = numel(s);
peak_mask = false(N, 1);
for n = 1:N
    left_val = -inf;
    right_val = -inf;
    if n > 1
        left_val = s(n - 1);
    end
    if n < N
        right_val = s(n + 1);
    end
    if s(n) >= left_val && s(n) >= right_val
        peak_mask(n) = true;
    end
end
idx = find(peak_mask);
if isempty(idx)
    [~, idx] = max(s);
end
[~, order] = sort(s(idx), 'descend');
idx = idx(order(1:min(numel(order), max_count)));
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
