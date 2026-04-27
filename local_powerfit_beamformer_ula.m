function [W, debug] = local_powerfit_beamformer_ula(M, beam_deg, design_grid_deg, d_eff, cfg)
% Localized power-domain beam synthesis.
% Compared with the previous global fit, this version emphasizes matching
% the mainlobe neighborhood and only softly suppresses the sidelobes.

if nargin < 5
    cfg = struct();
end

lambda_reg = get_opt(cfg, 'lambda_reg', 1e-2);
target_sigma_deg = get_opt(cfg, 'target_sigma_deg', 1.0);
target_floor = get_opt(cfg, 'target_floor', 0.10);
weight_sigma_deg = get_opt(cfg, 'weight_sigma_deg', 2.5);
far_weight = get_opt(cfg, 'far_weight', 0.12);
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
Wgt = zeros(P, G);
W = zeros(P, M);

for p = 1:P
    dist = design_grid_deg - beam_deg(p);
    main_env = exp(-0.5 * (dist / target_sigma_deg) .^ 2);
    target = target_floor + (1 - target_floor) * main_env;
    target = target / max(target);
    weight = far_weight + (1 - far_weight) * exp(-0.5 * (dist / weight_sigma_deg) .^ 2);

    T(p, :) = target;
    Wgt(p, :) = weight;

    amp_target = sqrt(target(:));
    phi = zeros(G, 1);
    for iter = 1:max_phase_iters
        b = amp_target .* exp(1i * phi);
        Aw = bsxfun(@times, A_H, sqrt(weight(:)));
        bw = sqrt(weight(:)) .* b;
        w = (Aw' * Aw + lambda_reg * eye(M)) \ (Aw' * bw);
        phi = angle(A_H * w);
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
debug.Wgt = Wgt;
debug.R = R;
debug.lambda_reg = lambda_reg;
debug.target_sigma_deg = target_sigma_deg;
debug.target_floor = target_floor;
debug.weight_sigma_deg = weight_sigma_deg;
debug.far_weight = far_weight;
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
