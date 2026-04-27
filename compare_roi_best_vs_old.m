function results = compare_roi_best_vs_old(mode_str)
if nargin < 1
    mode_str = 'full';
end

close all;
clc;
rng(20260415);

cfg = get_common_cfg(mode_str);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, ['compare_roi_best_vs_old_' timestamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

results.meta.cfg = cfg;
results.meta.output_dir = out_dir;
results.meta.timestamp = timestamp;

results.narrow = compare_one_case(get_narrow_case(), 'narrowband', cfg, ...
    struct('name', 'old-center_dense', 'shape', 'center_dense', 'roi_margin', []), ...
    struct('name', 'roi-center_dense-m4.25', 'shape', 'center_dense', 'roi_margin', 4.25));

results.wide = compare_one_case(get_wide_case(), 'wideband', cfg, ...
    struct('name', 'old-uniform', 'shape', 'uniform', 'roi_margin', []), ...
    struct('name', 'roi-uniform-m4.00', 'shape', 'uniform', 'roi_margin', 4.00));

write_summary(results, fullfile(out_dir, 'summary.txt'));
save_plots(results, out_dir);
save(fullfile(out_dir, 'results.mat'), 'results', '-v7.3');
fprintf('Saved focused ROI comparison to:\n%s\n', out_dir);
end

function case_result = compare_one_case(scene, modality, cfg, cand_a, cand_b)
candidate_list = {build_candidate(scene, modality, cfg, cand_a), build_candidate(scene, modality, cfg, cand_b)};
labels = cellfun(@(c) c.label, candidate_list, 'UniformOutput', false);
num_cand = numel(candidate_list);
num_snr = numel(cfg.snr_values);
mean_mse = zeros(num_cand, num_snr);

for c = 1:num_cand
    for s = 1:num_snr
        trial_mse = zeros(1, cfg.trials);
        for t = 1:cfg.trials
            data = generate_case_data(scene, modality, cfg.snr_values(s), candidate_list{c}.W, cfg, candidate_list{c}.freq_W_bands);
            if strcmp(modality, 'narrowband')
                [~, doa_hat] = zdoa_offgrid_refine(data.z, candidate_list{c}.W, cfg.scan_deg, scene.K, make_old_cfg(cfg));
            else
                [~, doa_hat] = zdoa_wideband_offgrid_refine(data.z_bands, data.D_bands, data.c_bands, data.W_bands, cfg.scan_deg, scene.K, make_old_cfg(cfg));
            end
            trial_mse(t) = best_match_mse(scene.true_doa_deg, doa_hat, cfg.fail_mse_deg2);
        end
        mean_mse(c, s) = mean(trial_mse);
    end
end

case_result.scene = scene;
case_result.modality = modality;
case_result.labels = labels;
case_result.snr_values = cfg.snr_values;
case_result.mean_mse_deg2 = mean_mse;
case_result.avg_mse_deg2 = mean(mean_mse, 2);
end

function candidate = build_candidate(scene, modality, cfg, spec)
if isempty(spec.roi_margin)
    beam_deg = design_beam_grid_local(cfg.P, cfg.scan_deg(1), cfg.scan_deg(end), spec.shape);
else
    roi_min = max(cfg.scan_deg(1), min(scene.true_doa_deg) - spec.roi_margin);
    roi_max = min(cfg.scan_deg(end), max(scene.true_doa_deg) + spec.roi_margin);
    beam_deg = design_beam_grid_local(cfg.P, roi_min, roi_max, spec.shape);
end
W = zeros(cfg.P, cfg.M);
for p = 1:cfg.P
    W(p, :) = steering_vector_local(beam_deg(p), cfg.M, cfg.d_over_lambda)';
end
W = normalize_rows(W);

candidate.label = spec.name;
candidate.W = W;
if strcmp(modality, 'wideband')
    candidate.freq_W_bands = repmat(W, 1, 1, numel(scene.freq_ratios));
else
    candidate.freq_W_bands = [];
end
end

function data = generate_case_data(scene, modality, snr_db, W, cfg, W_bands)
if nargin < 6
    W_bands = [];
end
if strcmp(modality, 'narrowband')
    X = generate_snapshots_band(cfg.M, cfg.L, scene.true_doa_deg, scene.source_power, snr_db, cfg.d_over_lambda);
    Y = W * X;
    data.z = mean(abs(Y) .^ 2, 2);
else
    num_bands = numel(scene.freq_ratios);
    z_bands = zeros(cfg.P, num_bands);
    D_bands = zeros(cfg.P, numel(cfg.scan_deg), num_bands);
    c_bands = zeros(cfg.P, num_bands);
    if isempty(W_bands)
        W_bands = repmat(W, 1, 1, num_bands);
    end
    for band_idx = 1:num_bands
        d_eff = cfg.d_over_lambda * scene.freq_ratios(band_idx);
        X_band = generate_snapshots_band(cfg.M, cfg.L, scene.true_doa_deg, scene.source_power, snr_db, d_eff);
        Y_band = W_bands(:, :, band_idx) * X_band;
        z_bands(:, band_idx) = mean(abs(Y_band) .^ 2, 2);
        D_bands(:, :, band_idx) = build_scan_dictionary(W_bands(:, :, band_idx), cfg.scan_deg, d_eff);
        c_bands(:, band_idx) = sum(abs(W_bands(:, :, band_idx)) .^ 2, 2);
    end
    data.z_bands = z_bands;
    data.D_bands = D_bands;
    data.c_bands = c_bands;
    data.W_bands = W_bands;
end
end

function D = build_scan_dictionary(W, scan_deg, d_eff)
D = zeros(size(W, 1), numel(scan_deg));
for idx = 1:numel(scan_deg)
    a = steering_vector_local(scan_deg(idx), size(W, 2), d_eff);
    D(:, idx) = abs(W * a) .^ 2;
end
end

function X = generate_snapshots_band(M, L, true_doa_deg, source_power, snr_db, d_eff)
K = numel(true_doa_deg);
A = zeros(M, K);
for k = 1:K
    A(:, k) = steering_vector_local(true_doa_deg(k), M, d_eff);
end
S = (randn(K, L) + 1i * randn(K, L)) / sqrt(2);
S = diag(sqrt(source_power(:))) * S;
signal = A * S;
signal_power = mean(abs(signal(:)) .^ 2);
noise_power = signal_power / (10 ^ (snr_db / 10));
noise = sqrt(noise_power / 2) * (randn(M, L) + 1i * randn(M, L));
X = signal + noise;
end

function beam_deg = design_beam_grid_local(P, scan_min, scan_max, design_name)
t = linspace(-1, 1, P);
switch design_name
    case 'uniform'
        u = t;
    case 'center_dense'
        u = sign(t) .* abs(t) .^ 1.8;
    case 'edge_dense'
        u = sign(t) .* abs(t) .^ 0.65;
    otherwise
        error('Unknown design_name');
end
beam_deg = (u + 1) / 2 * (scan_max - scan_min) + scan_min;
end

function W = normalize_rows(W)
row_norm = sqrt(sum(abs(W) .^ 2, 2));
W = bsxfun(@rdivide, W, max(row_norm, 1e-12));
end

function mse = best_match_mse(doa_true, doa_hat, fail_mse_deg2)
doa_true = sort(doa_true(:)).';
doa_hat = sort(doa_hat(:)).';
if numel(doa_true) ~= numel(doa_hat) || any(~isfinite(doa_hat))
    mse = fail_mse_deg2;
    return;
end
K = numel(doa_true);
perm_idx = perms(1:K);
best_val = inf;
for row_idx = 1:size(perm_idx, 1)
    err = doa_hat(perm_idx(row_idx, :)) - doa_true;
    val = mean(err .^ 2);
    if val < best_val
        best_val = val;
    end
end
mse = best_val;
end

function a = steering_vector_local(theta_deg, M, d_over_lambda)
n = (0:M - 1).';
a = exp(1i * 2 * pi * d_over_lambda * sind(theta_deg) * n) / sqrt(M);
end

function cfg = make_old_cfg(common)
cfg = struct();
cfg.L = common.L;
cfg.d_over_lambda = common.d_over_lambda;
cfg.fit_floor = 1e-8;
cfg.max_alt_iters = 12;
cfg.max_refine_passes = 2;
cfg.refine_radii_deg = [1.0, 0.5, 0.25, 0.12, 0.06];
cfg.refine_samples = 11;
cfg.refine_guard_deg = 0.6;
cfg.angle_min_deg = common.scan_deg(1);
cfg.angle_max_deg = common.scan_deg(end);
cfg.freq_ratios = common.freq_ratios;
cfg.max_support_evals = 50000;
end

function scene = get_narrow_case()
scene.name = 'NB-K3 close offgrid';
scene.K = 3;
scene.true_doa_deg = [-6.35, -2.05, 1.75];
scene.source_power = [1.0, 0.9, 0.7];
end

function scene = get_wide_case()
scene.name = 'WB-K2 strong-weak offgrid';
scene.K = 2;
scene.true_doa_deg = [-4.70, -0.70];
scene.source_power = [1.0, 0.33];
scene.freq_ratios = linspace(0.85, 1.15, 5);
end

function cfg = get_common_cfg(mode_str)
cfg.M = 20;
cfg.P = 20;
cfg.L = 40;
cfg.d_over_lambda = 0.5;
cfg.scan_deg = -18:0.5:18;
cfg.fail_mse_deg2 = 400;
cfg.freq_ratios = linspace(0.85, 1.15, 5);
switch lower(mode_str)
    case 'quick'
        cfg.snr_values = [-10, -6, -2, 2, 6, 10, 14];
        cfg.trials = 12;
    case 'full'
        cfg.snr_values = -10:2:16;
        cfg.trials = 32;
    otherwise
        error('mode_str must be quick or full.');
end
end

function write_summary(results, out_path)
fid = fopen(out_path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Focused ROI-best vs old-best comparison\n');
fprintf(fid, 'Timestamp: %s\n\n', results.meta.timestamp);
write_case(fid, results.narrow);
write_case(fid, results.wide);
end

function write_case(fid, case_result)
fprintf(fid, '[%s] %s\n', upper(case_result.modality), case_result.scene.name);
for i = 1:numel(case_result.labels)
    fprintf(fid, '  %s | avg MSE = %.4f deg^2\n', case_result.labels{i}, case_result.avg_mse_deg2(i));
end
fprintf(fid, '\n');
end

function save_plots(results, out_dir)
save_one_plot(results.narrow, fullfile(out_dir, 'narrow_compare.png'));
save_one_plot(results.wide, fullfile(out_dir, 'wide_compare.png'));
end

function save_one_plot(case_result, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 560]);
plot(case_result.snr_values, case_result.mean_mse_deg2(1, :), '-o', 'LineWidth', 2.2, 'MarkerSize', 7, 'Color', [0.85 0.33 0.10]); hold on;
plot(case_result.snr_values, case_result.mean_mse_deg2(2, :), '-s', 'LineWidth', 2.2, 'MarkerSize', 7, 'Color', [0 0.45 0.74]);
grid on;
xlabel('SNR (dB)');
ylabel('MSE (deg^2)');
title(sprintf('%s: ROI best vs old best', case_result.scene.name), 'Interpreter', 'none');
legend(case_result.labels, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 13);
exportgraphics(fig, out_path, 'Resolution', 220);
close(fig);
end
