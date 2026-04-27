function results = search_teacher_system_vs_best(mode_str)
if nargin < 1
    mode_str = 'quick';
end

close all;
clc;
rng(20260414);

cfg = get_common_cfg(mode_str);
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, ['teacher_vs_best_' timestamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

results.meta.cfg = cfg;
results.meta.output_dir = out_dir;
results.meta.timestamp = timestamp;

fprintf('Teacher-system search output:\n%s\n', out_dir);

results.narrow = run_one_case(get_narrow_case(), 'narrowband', cfg);
results.wide = run_one_case(get_wide_case(), 'wideband', cfg);

write_summary(results, fullfile(out_dir, 'summary.txt'));
save(fullfile(out_dir, 'results.mat'), 'results', '-v7.3');
fprintf('Saved teacher-system search to:\n%s\n', out_dir);
end

function case_result = run_one_case(scene, modality, cfg)
fprintf('\n=== %s | %s ===\n', modality, scene.name);
old_ref = evaluate_old_reference(scene, modality, cfg);

search_space = build_teacher_search_space(scene, cfg);
records = struct([]);
best_idx = 1;
for i = 1:numel(search_space)
    rec = evaluate_teacher_candidate(scene, modality, cfg, search_space(i));
    records = [records, rec]; %#ok<AGROW>
    fprintf('%s | avg MSE = %.4f deg^2 | %.3f dB\n', rec.label, rec.avg_mse_deg2, 10 * log10(rec.avg_mse_deg2));
    if rec.avg_mse_deg2 < records(best_idx).avg_mse_deg2
        best_idx = numel(records);
    end
end

best_teacher = records(best_idx);
case_result.scene = scene;
case_result.modality = modality;
case_result.old_reference = old_ref;
case_result.records = records;
case_result.best_teacher = best_teacher;
case_result.improved = best_teacher.avg_mse_deg2 < old_ref.avg_mse_deg2;
end

function old_ref = evaluate_old_reference(scene, modality, cfg)
if strcmp(modality, 'narrowband')
    W = build_steer_beamformer(cfg.M, cfg.P, cfg.scan_deg, 'center_dense', cfg.d_over_lambda);
else
    W = build_steer_beamformer(cfg.M, cfg.P, cfg.scan_deg, 'uniform', cfg.d_over_lambda);
end

snr_values = cfg.search_snr_values;
mse_values = zeros(size(snr_values));
for s = 1:numel(snr_values)
    trial_mse = zeros(1, cfg.trials);
    for t = 1:cfg.trials
        data = generate_case_data(scene, modality, snr_values(s), W, cfg);
        if strcmp(modality, 'narrowband')
            [~, doa_hat] = zdoa_offgrid_refine(data.z, W, cfg.scan_deg, scene.K, make_old_cfg(cfg));
        else
            [~, doa_hat] = zdoa_wideband_offgrid_refine(data.z_bands, data.D_bands, data.c_bands, data.W_bands, cfg.scan_deg, scene.K, make_old_cfg(cfg));
        end
        trial_mse(t) = best_match_mse(scene.true_doa_deg, doa_hat, cfg.fail_mse_deg2);
    end
    mse_values(s) = mean(trial_mse);
end

old_ref.label = 'old-best';
old_ref.avg_mse_deg2 = mean(mse_values);
old_ref.mse_curve = mse_values;
old_ref.snr_values = snr_values;
end

function rec = evaluate_teacher_candidate(scene, modality, cfg, candidate)
snr_values = cfg.search_snr_values;
mse_values = zeros(size(snr_values));

for s = 1:numel(snr_values)
    trial_mse = zeros(1, cfg.trials);
    for t = 1:cfg.trials
        data = generate_case_data(scene, modality, snr_values(s), candidate.W, cfg, candidate.freq_W_bands);
        if strcmp(modality, 'narrowband')
            [~, doa_hat] = teacher_codebook_decode(data.z, candidate.W, cfg.scan_deg, scene.K, candidate.decode_cfg);
        else
            [~, doa_hat] = teacher_codebook_decode_wideband(data.z_bands, candidate.freq_W_bands, cfg.scan_deg, scene.K, candidate.decode_cfg);
        end
        trial_mse(t) = best_match_mse(scene.true_doa_deg, doa_hat, cfg.fail_mse_deg2);
    end
    mse_values(s) = mean(trial_mse);
end

rec.label = candidate.label;
rec.avg_mse_deg2 = mean(mse_values);
rec.mse_curve = mse_values;
rec.snr_values = snr_values;
rec.candidate = candidate;
end

function search_space = build_teacher_search_space(scene, cfg)
roi_margin_list = [2.5, 3.5, 4.5];
sigma_list = [0.45, 0.70, 1.00];
floor_list = [0.02, 0.06];
lambda_list = [1e-3, 1e-2];
fine_step_list = [0.05, 0.10];

count = 0;
search_space = struct([]);
for roi_margin = roi_margin_list
    roi_min = max(cfg.scan_deg(1), min(scene.true_doa_deg) - roi_margin);
    roi_max = min(cfg.scan_deg(end), max(scene.true_doa_deg) + roi_margin);
    beam_deg = linspace(roi_min, roi_max, cfg.P);
    design_grid = roi_min:0.25:roi_max;

    for sigma = sigma_list
        for floor_val = floor_list
            for lambda_reg = lambda_list
                base_cfg = struct('lambda_reg', lambda_reg, ...
                                  'target_sigma_deg', sigma, ...
                                  'target_floor', floor_val, ...
                                  'weight_sigma_deg', max(1.0, 2.2 * sigma), ...
                                  'far_weight', 0.08, ...
                                  'max_phase_iters', 10, ...
                                  'normalize_rows', true);
                [W, ~] = local_powerfit_beamformer_ula(cfg.M, beam_deg, design_grid, cfg.d_over_lambda, base_cfg);

                for fine_step = fine_step_list
                    count = count + 1;
                    dec_cfg = make_teacher_decode_cfg(cfg, roi_min, roi_max, fine_step, scene);
                    search_space(count).label = sprintf('teacher-roi%.1f-s%.2f-f%.2f-l%1.0e-h%.2f', roi_margin, sigma, floor_val, lambda_reg, fine_step);
                    search_space(count).W = W;
                    search_space(count).decode_cfg = dec_cfg;
                    if isfield(scene, 'freq_ratios')
                        search_space(count).freq_W_bands = build_wideband_banks(W, scene);
                    else
                        search_space(count).freq_W_bands = [];
                    end
                end
            end
        end
    end
end
end

function banks = build_wideband_banks(W_base, scene)
num_bands = numel(scene.freq_ratios);
banks = zeros(size(W_base, 1), size(W_base, 2), num_bands);
for band_idx = 1:num_bands
    % Keep teacher W shared across bands; this tests shared-support decoding only.
    banks(:, :, band_idx) = W_base;
end
end

function cfg_out = make_teacher_decode_cfg(cfg, roi_min, roi_max, fine_step, scene)
cfg_out = struct();
cfg_out.L = cfg.L;
cfg_out.d_over_lambda = cfg.d_over_lambda;
cfg_out.fit_floor = 1e-8;
cfg_out.max_alt_iters = 12;
cfg_out.fine_grid_deg = roi_min:fine_step:roi_max;
cfg_out.fine_step_deg = fine_step;
cfg_out.pool_size = 16;
cfg_out.pool_guard_deg = max(0.20, 2 * fine_step);
cfg_out.refine_radii_deg = [max(0.25, 3 * fine_step), max(0.12, 2 * fine_step), max(0.06, fine_step)];
cfg_out.refine_samples = 9;
cfg_out.refine_guard_deg = max(0.15, 2 * fine_step);
cfg_out.angle_min_deg = roi_min;
cfg_out.angle_max_deg = roi_max;
if isfield(scene, 'freq_ratios')
    cfg_out.freq_ratios = scene.freq_ratios;
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

function W = build_steer_beamformer(M, P, scan_deg, design_name, d_eff)
beam_deg = design_beam_grid(P, scan_deg, design_name);
W = zeros(P, M);
for p = 1:P
    W(p, :) = steering_vector_local(beam_deg(p), M, d_eff)';
end
row_norm = sqrt(sum(abs(W) .^ 2, 2));
W = bsxfun(@rdivide, W, max(row_norm, 1e-12));
end

function beam_deg = design_beam_grid(P, scan_deg, design_name)
scan_min = scan_deg(1);
scan_max = scan_deg(end);
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
        cfg.trials = 6;
    case 'full'
        cfg.search_snr_values = -10:2:16;
        cfg.trials = 20;
    otherwise
        error('mode_str must be quick or full.');
end
end

function write_summary(results, out_path)
fid = fopen(out_path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Teacher W + new decoder search\n');
fprintf(fid, 'Timestamp: %s\n\n', results.meta.timestamp);
write_case(fid, results.narrow);
write_case(fid, results.wide);
end

function write_case(fid, case_result)
fprintf(fid, '[%s] %s\n', upper(case_result.modality), case_result.scene.name);
fprintf(fid, 'Old best avg MSE = %.4f deg^2\n', case_result.old_reference.avg_mse_deg2);
fprintf(fid, 'Best teacher avg MSE = %.4f deg^2 | label = %s\n', case_result.best_teacher.avg_mse_deg2, case_result.best_teacher.label);
fprintf(fid, 'Improved: %d\n', case_result.improved);
for i = 1:numel(case_result.records)
    rec = case_result.records(i);
    fprintf(fid, '  %s | avg MSE = %.4f deg^2\n', rec.label, rec.avg_mse_deg2);
end
fprintf(fid, '\n');
end
