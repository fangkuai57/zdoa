function [W, debug] = powerfit_beamformer_ula(M, beam_deg, design_grid_deg, d_eff, cfg)
% Power-domain beam synthesis.
% We design desired power templates, then fit complex amplitudes iteratively:
%   |A^H w|^2 ~= T
% by alternating between phase update and regularized LS on amplitude.

if nargin < 5
    cfg = struct();
end

lambda_reg = get_opt(cfg, 'lambda_reg', 5e-3);
target_sigma_deg = get_opt(cfg, 'target_sigma_deg', 1.1);
target_floor = get_opt(cfg, 'target_floor', 0.08);
target_clip = get_opt(cfg, 'target_clip', 0.98);
max_phase_iters = get_opt(cfg, 'max_phase_iters', 8);
normalize_rows = get_opt(cfg, 'normalize_rows', true);

G = numel(design_grid_deg);
P = numel(beam_deg);

A = zeros(M, G);
for g = 1:G
    A(:, g) = steering_vector_local(design_grid_deg(g), M, d_eff);
end
A_H = A';

T = zeros(P, G);
for p = 1:P
    dist = design_grid_deg - beam_deg(p);
    target = target_floor + (1 - target_floor) * exp(-0.5 * (dist / target_sigma_deg) .^ 2);
    target = min(target, target_clip);
    target = target / max(target);
    T(p, :) = target;
end

W = zeros(P, M);
amp_target = sqrt(T);
for p = 1:P
    phi = zeros(G, 1);
    w = zeros(M, 1);
    for iter = 1:max_phase_iters
        b = amp_target(p, :).'.* exp(1i * phi);
        w = (A * A_H + lambda_reg * eye(M)) \ (A * b);
        response = A_H * w;
        phi = angle(response);
    end
    W(p, :) = w.';
end

if normalize_rows
    row_norm = sqrt(sum(abs(W) .^ 2, 2));
    row_norm = max(row_norm, 1e-12);
    W = bsxfun(@rdivide, W, row_norm);
end

R = abs(W * A) .^ 2;
debug.A = A;
debug.T = T;
debug.R = R;
debug.lambda_reg = lambda_reg;
debug.target_sigma_deg = target_sigma_deg;
debug.target_floor = target_floor;
debug.beam_deg = beam_deg(:).';
debug.design_grid_deg = design_grid_deg(:).';
end

function value = get_opt(cfg, field_name, default_value)
if isfield(cfg, field_name)
    value = cfg.(field_name);
else
    value = default_value;
end
end

function a = steering_vector_local(theta_deg, M, d_eff)
n = (0:M - 1).';
a = exp(1i * 2 * pi * d_eff * sind(theta_deg) * n) / sqrt(M);
end
