function [W, debug] = impulse_fit_beamformer_ula(M, beam_deg, design_grid_deg, d_eff, cfg)
% Regularized beam synthesis from desired angular response.
% Solve W * A_design ~= B_des in least-squares sense.

if nargin < 5
    cfg = struct();
end

lambda_reg = get_opt(cfg, 'lambda_reg', 1e-2);
target_sigma_deg = get_opt(cfg, 'target_sigma_deg', 0.55);
target_floor = get_opt(cfg, 'target_floor', 0.0);
normalize_rows = get_opt(cfg, 'normalize_rows', true);

G = numel(design_grid_deg);
P = numel(beam_deg);

A = zeros(M, G);
for g = 1:G
    A(:, g) = steering_vector_local(design_grid_deg(g), M, d_eff);
end

B_des = zeros(P, G);
for p = 1:P
    dist = design_grid_deg - beam_deg(p);
    B_des(p, :) = exp(-0.5 * (dist / target_sigma_deg) .^ 2);
    if target_floor > 0
        B_des(p, B_des(p, :) < target_floor) = 0;
    end
    row_peak = max(B_des(p, :));
    if row_peak > 0
        B_des(p, :) = B_des(p, :) / row_peak;
    end
end

W = B_des * A' / (A * A' + lambda_reg * eye(M));

if normalize_rows
    row_norm = sqrt(sum(abs(W) .^ 2, 2));
    row_norm = max(row_norm, 1e-12);
    W = bsxfun(@rdivide, W, row_norm);
end

debug.A = A;
debug.B_des = B_des;
debug.lambda_reg = lambda_reg;
debug.target_sigma_deg = target_sigma_deg;
debug.design_grid_deg = design_grid_deg(:).';
debug.beam_deg = beam_deg(:).';
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
