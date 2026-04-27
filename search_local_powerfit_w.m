function results = search_local_powerfit_w(mode_str)
if nargin < 1
    mode_str = 'quick';
end

close all;
clc;
rng(20260414);

cfg = get_cfg(mode_str);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, ['search_local_powerfit_w_' timestamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

results.meta = cfg;
results.meta.output_dir = out_dir;
results.meta.timestamp = timestamp;

results.narrow = run_search(make_case('narrowband', cfg), cfg);
results.wide = run_search(make_case('wideband', cfg), cfg);

write_summary(results, fullfile(out_dir, 'summary.txt'));
save(fullfile(out_dir, 'results.mat'), 'results', '-v7.3');
fprintf('Saved local-powerfit W search to:\n%s\n', out_dir);
end

function cfg = get_cfg(mode_str)
cfg.M = 20;
cfg.P = 20;
cfg.L = 40;
cfg.scan_deg = -18:0.5:18;
cfg.snr_values = [-10, -2, 6, 14];
cfg.d_over_lambda = 0.5;
cfg.max_support_evals = 50000;
cfg.max_alt_iters = 12;
cfg.fit_floor = 1e-8;
cfg.refine_radii_deg = [1.0, 0.5, 0.25, 0.12, 0.06];
cfg.refine_samples = 11;
cfg.refine_guard_deg = 0.6;
cfg.freq_ratios = [0.85 0.925 1.0 1.075 1.15];
switch lower(mode_str)
    case 'quick'
        cfg.trials = 2;
    case 'full'
        cfg.trials = 6;
    otherwise
        error('mode_str must be quick or full.');
end
end

function case_cfg = make_case(modality, cfg)
case_cfg = cfg;
case_cfg.modality = modality;
if strcmp(modality, 'narrowband')
    case_cfg.name = 'NB-K3 close offgrid';
    case_cfg.true_doa_deg = [-6.35, -2.05, 1.75];
    case_cfg.source_power = [1.0, 0.9, 0.7];
    case_cfg.K = 3;
else
    case_cfg.name = 'WB-K2 strong-weak offgrid';
    case_cfg.true_doa_deg = [-4.70, -0.70];
    case_cfg.source_power = [1.0, 0.33];
    case_cfg.K = 2;
end
end

function result = run_search(case_cfg, cfg)
cands = build_candidates();
result.case_cfg = case_cfg;
result.candidates = cands;
result.avg_mse_deg2 = zeros(numel(cands), 1);
result.avg_mse_db = zeros(numel(cands), 1);
result.coherence = zeros(numel(cands), 1);

for idx = 1:numel(cands)
    frontend = build_frontend(case_cfg, cands{idx});
    result.coherence(idx) = dictionary_coherence(frontend.D);
    mse_total = 0;
    count = 0;
    for s = 1:numel(cfg.snr_values)
        for trial = 1:cfg.trials
            data = generate_data(case_cfg, frontend, cfg.snr_values(s));
            doa_hat = estimate_one(data, frontend, case_cfg);
            mse_total = mse_total + best_match_mse(sort(doa_hat(:)).', sort(case_cfg.true_doa_deg(:)).');
            count = count + 1;
        end
    end
    result.avg_mse_deg2(idx) = mse_total / count;
    result.avg_mse_db(idx) = 10 * log10(max(result.avg_mse_deg2(idx), 1e-4));
    fprintf('%s | %s | avg MSE = %.4f deg^2 | %.3f dB | coh = %.4f\n', ...
        case_cfg.modality, cands{idx}.label, result.avg_mse_deg2(idx), result.avg_mse_db(idx), result.coherence(idx));
end

[~, best_idx] = min(result.avg_mse_deg2);
result.best_idx = best_idx;
result.best_candidate = cands{best_idx};
end

function cands = build_candidates()
shapes = {'uniform', 'center_dense', 'edge_dense'};
cands = {};
for s = 1:numel(shapes)
    cands{end + 1} = struct('kind', 'steer', 'beam_shape', shapes{s}, 'label', ['steer-' shapes{s}]); %#ok<AGROW>
end
sigmas = [0.8, 1.2];
floors = [0.05, 0.12];
lambdas = [1e-3, 1e-2];
for s = 1:numel(shapes)
    for g = 1:numel(sigmas)
        for f = 1:numel(floors)
            for l = 1:numel(lambdas)
                cands{end + 1} = struct( ...
                    'kind', 'local_powerfit', ...
                    'beam_shape', shapes{s}, ...
                    'target_sigma_deg', sigmas(g), ...
                    'target_floor', floors(f), ...
                    'lambda_reg', lambdas(l), ...
                    'weight_sigma_deg', 2.5, ...
                    'far_weight', 0.12, ...
                    'label', sprintf('localpf-%s-s%.1f-f%.2f-l%.0e', shapes{s}, sigmas(g), floors(f), lambdas(l))); %#ok<AGROW>
            end
        end
    end
end
end

function frontend = build_frontend(case_cfg, cand)
beam_deg = design_beam_grid(case_cfg.P, case_cfg.scan_deg, cand.beam_shape);
num_bands = numel(case_cfg.freq_ratios);
frontend.W_bands = zeros(case_cfg.P, case_cfg.M, num_bands);
frontend.D_bands = zeros(case_cfg.P, numel(case_cfg.scan_deg), num_bands);
frontend.c_bands = zeros(case_cfg.P, num_bands);

for b = 1:num_bands
    d_eff = case_cfg.d_over_lambda * case_cfg.freq_ratios(b);
    if strcmp(cand.kind, 'steer')
        W = build_steered_beamformer_local(case_cfg.M, beam_deg, d_eff);
    else
        [W, ~] = local_powerfit_beamformer_ula(case_cfg.M, beam_deg, case_cfg.scan_deg, d_eff, struct( ...
            'lambda_reg', cand.lambda_reg, ...
            'target_sigma_deg', cand.target_sigma_deg, ...
            'target_floor', cand.target_floor, ...
            'weight_sigma_deg', cand.weight_sigma_deg, ...
            'far_weight', cand.far_weight));
    end
    [D, c] = build_beam_power_dictionary_local(W, case_cfg.scan_deg, d_eff);
    frontend.W_bands(:, :, b) = W;
    frontend.D_bands(:, :, b) = D;
    frontend.c_bands(:, b) = c;
end
frontend.W = frontend.W_bands(:, :, 1);
frontend.D = frontend.D_bands(:, :, 1);
frontend.c = frontend.c_bands(:, 1);
end

function data = generate_data(case_cfg, frontend, snr_db)
if strcmp(case_cfg.modality, 'narrowband')
    X = generate_snapshots_band(case_cfg.M, case_cfg.L, case_cfg.true_doa_deg, case_cfg.source_power, snr_db, case_cfg.d_over_lambda);
    Y = frontend.W * X;
    data.z = mean(abs(Y) .^ 2, 2);
else
    num_bands = numel(case_cfg.freq_ratios);
    z_bands = zeros(case_cfg.P, num_bands);
    for b = 1:num_bands
        d_eff = case_cfg.d_over_lambda * case_cfg.freq_ratios(b);
        X = generate_snapshots_band(case_cfg.M, case_cfg.L, case_cfg.true_doa_deg, case_cfg.source_power, snr_db, d_eff);
        Y = frontend.W_bands(:, :, b) * X;
        z_bands(:, b) = mean(abs(Y) .^ 2, 2);
    end
    data.z_bands = z_bands;
end
end

function doa_hat = estimate_one(data, frontend, case_cfg)
algo_cfg = struct( ...
    'L', case_cfg.L, ...
    'd_over_lambda', case_cfg.d_over_lambda, ...
    'freq_ratios', case_cfg.freq_ratios, ...
    'max_support_evals', case_cfg.max_support_evals, ...
    'max_alt_iters', case_cfg.max_alt_iters, ...
    'fit_floor', case_cfg.fit_floor, ...
    'refine_radii_deg', case_cfg.refine_radii_deg, ...
    'refine_samples', case_cfg.refine_samples, ...
    'refine_guard_deg', case_cfg.refine_guard_deg, ...
    'angle_min_deg', case_cfg.scan_deg(1), ...
    'angle_max_deg', case_cfg.scan_deg(end));
if strcmp(case_cfg.modality, 'narrowband')
    [~, doa_hat] = zdoa_offgrid_refine(data.z, frontend.W, case_cfg.scan_deg, case_cfg.K, algo_cfg);
else
    [~, doa_hat] = zdoa_wideband_offgrid_refine(data.z_bands, frontend.D_bands, frontend.c_bands, frontend.W_bands, case_cfg.scan_deg, case_cfg.K, algo_cfg);
end
end

function write_summary(results, filename)
fid = fopen(filename, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Local-powerfit W search\n');
fprintf(fid, 'Timestamp: %s\n\n', results.meta.timestamp);
write_one(fid, 'NARROW', results.narrow);
fprintf(fid, '\n');
write_one(fid, 'WIDE', results.wide);
end

function write_one(fid, title_str, result)
fprintf(fid, '[%s] %s\n', title_str, result.case_cfg.name);
fprintf(fid, 'Best candidate: %s\n', result.best_candidate.label);
for idx = 1:numel(result.candidates)
    fprintf(fid, '  %-34s avg MSE = %.4f deg^2 | %.3f dB | coh = %.4f\n', ...
        result.candidates{idx}.label, ...
        result.avg_mse_deg2(idx), ...
        result.avg_mse_db(idx), ...
        result.coherence(idx));
end
end

function coh = dictionary_coherence(D)
D = D ./ max(vecnorm(D, 2, 1), 1e-12);
G = abs(D' * D);
G(1:size(G, 1)+1:end) = 0;
coh = max(G(:));
end

function beam_deg = design_beam_grid(P, scan_deg, design_name)
scan_min = scan_deg(1);
scan_max = scan_deg(end);
t = linspace(-1, 1, P);
switch design_name
    case 'uniform'
        beam_deg = linspace(scan_min, scan_max, P);
    case 'center_dense'
        beam_deg = max(abs([scan_min, scan_max])) * sign(t) .* abs(t) .^ 1.5;
    case 'edge_dense'
        beam_deg = max(abs([scan_min, scan_max])) * sin((pi / 2) * t);
    otherwise
        error('Unknown design_name %s', design_name);
end
beam_deg = sort(max(min(beam_deg, scan_max), scan_min));
end

function [D, c] = build_beam_power_dictionary_local(W, scan_deg, d_eff)
P = size(W, 1);
G = numel(scan_deg);
M = size(W, 2);
D = zeros(P, G);
for g = 1:G
    a = steering_vector_local(scan_deg(g), M, d_eff);
    D(:, g) = abs(W * a) .^ 2;
end
c = sum(abs(W) .^ 2, 2);
end

function W = build_steered_beamformer_local(M, beam_deg, d_eff)
W = zeros(numel(beam_deg), M);
for p = 1:numel(beam_deg)
    W(p, :) = steering_vector_local(beam_deg(p), M, d_eff)';
end
end

function X = generate_snapshots_band(M, L, true_doa_deg, source_power, snr_db, d_eff)
K = numel(true_doa_deg);
A = zeros(M, K);
for k = 1:K
    A(:, k) = steering_vector_local(true_doa_deg(k), M, d_eff);
end
S = (randn(K, L) + 1i * randn(K, L)) / sqrt(2);
for k = 1:K
    S(k, :) = sqrt(source_power(k)) * S(k, :);
end
signal_var = mean(abs(A * S) .^ 2, 'all');
noise_var = signal_var / (10 ^ (snr_db / 10));
N = sqrt(noise_var / 2) * (randn(M, L) + 1i * randn(M, L));
X = A * S + N;
end

function a = steering_vector_local(theta_deg, M, d_eff)
n = (0:M - 1).';
a = exp(1i * 2 * pi * d_eff * sind(theta_deg) * n) / sqrt(M);
end

function value = best_match_mse(doa_hat, doa_true)
perm_idx = perms(1:numel(doa_true));
best_val = inf;
for row = 1:size(perm_idx, 1)
    err = doa_hat - doa_true(perm_idx(row, :));
    best_val = min(best_val, mean(err .^ 2));
end
value = best_val;
end
