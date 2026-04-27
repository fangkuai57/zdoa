function [W, debug] = teacher_inverse_beamformer_ula(M, design_points_deg, d_eff, cfg)
% Teacher-style inverse beamformer:
% choose desired discrete angular responses first, then solve W*A ~= I.

if nargin < 4
    cfg = struct();
end

lambda_reg = get_opt(cfg, 'lambda_reg', 1e-2);
normalize_rows = get_opt(cfg, 'normalize_rows', true);
embed_scale = get_opt(cfg, 'embed_scale', 1.0);

G = numel(design_points_deg);
A = zeros(M, G);
for g = 1:G
    A(:, g) = steering_vector_local(design_points_deg(g), M, d_eff);
end

B = embed_scale * eye(G);
W = B * A' / (A * A' + lambda_reg * eye(M));

if normalize_rows
    row_norm = sqrt(sum(abs(W) .^ 2, 2));
    row_norm = max(row_norm, 1e-12);
    W = bsxfun(@rdivide, W, row_norm);
end

R = abs(W * A) .^ 2;
debug.design_points_deg = design_points_deg(:).';
debug.response = R;
debug.lambda_reg = lambda_reg;
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
