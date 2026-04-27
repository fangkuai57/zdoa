function results = search_roi_margin_refine_vs_best(mode_str)
if nargin < 1
    mode_str = 'quick';
end

close all;
clc;
rng(20260415);

cfg = get_common_cfg(mode_str);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, ['roi_margin_refine_vs_best_' timestamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

results.meta.cfg = cfg;
results.meta.output_dir = out_dir;
results.meta.timestamp = timestamp;

results.narrow = run_case(get_narrow_case(), 'narrowband', cfg, 'center_dense', [3.0, 3.25, 3.5, 3.75, 4.0, 4.25]);
results.wide = run_case(get_wide_case(), 'wideband', cfg, 'uniform', [3.0, 3.25, 3.5, 3.75, 4.0, 4.25]);

write_summary(results, fullfile(out_dir, 'summary.txt'));
save(fullfile(out_dir, 'results.mat'), 'results', '-v7.3');
fprintf('Saved ROI margin refine search to:\n%s\n', out_dir);
end

function case_result = run_case(scene, modality, cfg, shape_name, margin_list)
fprintf('\n=== %s | %s ===\n', modality, scene.name);
old_ref = evaluate_candidate(scene, modality, cfg, build_candidate(scene, modality, cfg, shape_name, []));

records = {};
for i = 1:numel(margin_list)
    records{i} = build_candidate(scene, modality, cfg, shape_name, margin_list(i)); %#ok<AGROW>
    eval_out = evaluate_candidate(scene, modality, cfg, records{i});
    records{i}.avg_mse_deg2 = eval_out.avg_mse_deg2;
    records{i}.mse_curve = eval_out.mse_curve;
    fprintf('%s | avg MSE = %.4f deg^2 | %.3f dB\n', records{i}.label, records{i}.avg_mse_deg2, 10 * log10(records{i}.avg_mse_deg2));
end

[~, best_idx] = min(cellfun(@(r) r.avg_mse_deg2, records));
case_result.scene = scene;
case_result.modality = modality;
case_result.old_reference = old_ref;
case_result.records = records;
case_result.best_teacher = records{best_idx};
case_result.improved = records{best_idx}.avg_mse_deg2 < old_ref.avg_mse_deg2;
end

function candidate = build_candidate(scene, modality, cfg, shape_name, roi_margin)
if isempty(roi_margin)
    beam_deg = design_beam_grid_local(cfg.P, cfg.scan_deg(1), cfg.scan_deg(end), shape_name);
    label = ['old-' shape_name];
else
    roi_min = max(cfg.scan_deg(1), min(scene.true_doa_deg) - roi_margin);
    roi_max = min(cfg.scan_deg(end), max(scene.true_doa_deg) + roi_margin);
    beam_deg = design_beam_grid_local(cfg.P, roi_min, roi_max, shape_name);
    label = sprintf('roi-%s-m%.2f', shape_name, roi_margin);
end

W = zeros(cfg.P, cfg.M);
for p = 1:cfg.P
    W(p, :) = steering_vector_local(beam_deg(p), cfg.M, cfg.d_over_lambda)';
end
W = normalize_rows(W);

candidate.label = label;
candidate.W = W;
if strcmp(modality, 'wideband')
    candidate.freq_W_bands = repmat(W, 1, 1, numel(scene.freq_ratios));
else
    candidate.freq_W_bands = [];
end
end

function out = evaluate_candidate(scene, modality, cfg, candidate)
snr_values = cfg.search_snr_values;
mse_values = zeros(size(snr_values));
for s = 1:numel(snr_values)
    trial_mse = zeros(1, cfg.trials);
    for t = 1:cfg.trials
        data = generate_case_data(scene, modality, snr_values(s), candidate.W, cfg, candidate.freq_W_bands);
        if strcmp(modality, 'narrowband')
            [~, doa_hat] = zdoa_offgrid_refine(data.z, candidate.W, cfg.scan_deg, scene.K, make_old_cfg(cfg));
        else
            [~, doa_hat] = zdoa_wideband_offgrid_refine(data.z_bands, data.D_bands, data.c_bands, data.W_bands, cfg.scan_deg, scene.K, make_old_cfg(cfg));
        end
        trial_mse(t) = best_match_mse(scene.true_doa_deg, doa_hat, cfg.fail_mse_deg2);
    end
    mse_values(s) = mean(trial_mse);
end
out.avg_mse_deg2 = mean(mse_values);
out.mse_curve = mse_values;
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
        cfg.search_snr_values = [-10, -6, -2, 2, 6, 10, 14];
        cfg.trials = 12;
    case 'full'
        cfg.search_snr_values = -10:2:16;
        cfg.trials = 32;
    otherwise
        error('mode_str must be quick or full.');
end
end

function write_summary(results, out_path)
fid = fopen(out_path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'ROI margin refine search\n');
fprintf(fid, 'Timestamp: %s\n\n', results.meta.timestamp);
write_case(fid, results.narrow);
write_case(fid, results.wide);
end

function write_case(fid, case_result)
fprintf(fid, '[%s] %s\n', upper(case_result.modality), case_result.scene.name);
fprintf(fid, 'Old best avg MSE = %.4f deg^2\n', case_result.old_reference.avg_mse_deg2);
fprintf(fid, 'Best refined avg MSE = %.4f deg^2 | label = %s\n', case_result.best_teacher.avg_mse_deg2, case_result.best_teacher.label);
fprintf(fid, 'Improved: %d\n', case_result.improved);
for i = 1:numel(case_result.records)
    rec = case_result.records{i};
    fprintf(fid, '  %s | avg MSE = %.4f deg^2\n', rec.label, rec.avg_mse_deg2);
end
fprintf(fid, '\n');
end
