function results = search_wide_roi_margin_medium()
close all;
clc;
rng(20260416);

cfg.M = 20;
cfg.P = 20;
cfg.L = 40;
cfg.d_over_lambda = 0.5;
cfg.scan_deg = -18:0.5:18;
cfg.fail_mse_deg2 = 400;
cfg.freq_ratios = linspace(0.85, 1.15, 5);
cfg.snr_values = -10:2:16;
cfg.trials = 20;

scene.name = 'WB-K2 strong-weak offgrid';
scene.K = 2;
scene.true_doa_deg = [-4.70, -0.70];
scene.source_power = [1.0, 0.33];
scene.freq_ratios = cfg.freq_ratios;

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
out_dir = fullfile(pwd, ['wide_roi_margin_medium_' timestamp]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

margins = [3.6, 3.8, 4.0, 4.2, 4.4, 4.6];
labels = ['old', compose('roi-m%.1f', margins)];
mean_mse = zeros(numel(labels), numel(cfg.snr_values));

for c = 1:numel(labels)
    if c == 1
        beam_deg = design_beam_grid_local(cfg.P, cfg.scan_deg(1), cfg.scan_deg(end), 'uniform');
    else
        margin = margins(c - 1);
        roi_min = max(cfg.scan_deg(1), min(scene.true_doa_deg) - margin);
        roi_max = min(cfg.scan_deg(end), max(scene.true_doa_deg) + margin);
        beam_deg = design_beam_grid_local(cfg.P, roi_min, roi_max, 'uniform');
    end

    W = build_W_from_beams(beam_deg, cfg.M, cfg.d_over_lambda);
    W_bands = repmat(W, 1, 1, numel(scene.freq_ratios));

    for s = 1:numel(cfg.snr_values)
        trial_mse = zeros(1, cfg.trials);
        for t = 1:cfg.trials
            data = generate_case_data(scene, cfg, cfg.snr_values(s), W_bands);
            [~, doa_hat] = zdoa_wideband_offgrid_refine(data.z_bands, data.D_bands, data.c_bands, data.W_bands, cfg.scan_deg, scene.K, make_old_cfg(cfg));
            trial_mse(t) = best_match_mse(scene.true_doa_deg, doa_hat, cfg.fail_mse_deg2);
        end
        mean_mse(c, s) = mean(trial_mse);
    end
    fprintf('%s | avg MSE = %.4f deg^2 | %.3f dB\n', labels{c}, mean(mean_mse(c, :)), 10 * log10(mean(mean_mse(c, :))));
end

[best_val, best_idx] = min(mean(mean_mse, 2));
results.labels = labels;
results.mean_mse_deg2 = mean_mse;
results.avg_mse_deg2 = mean(mean_mse, 2);
results.best_label = labels{best_idx};
results.best_mse_deg2 = best_val;
results.output_dir = out_dir;

fid = fopen(fullfile(out_dir, 'summary.txt'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Wide ROI margin medium search\n');
for c = 1:numel(labels)
    fprintf(fid, '%s | avg MSE = %.4f deg^2\n', labels{c}, results.avg_mse_deg2(c));
end
fprintf(fid, 'Best: %s | %.4f deg^2\n', results.best_label, results.best_mse_deg2);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 920, 560]);
colors = lines(numel(labels));
for c = 1:numel(labels)
    plot(cfg.snr_values, mean_mse(c, :), '-o', 'LineWidth', 1.8, 'MarkerSize', 6, 'Color', colors(c, :)); hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('MSE (deg^2)');
title('Wide ROI margin medium search');
legend(labels, 'Location', 'northeast');
exportgraphics(fig, fullfile(out_dir, 'wide_roi_margin_medium.png'), 'Resolution', 220);
close(fig);

save(fullfile(out_dir, 'results.mat'), 'results', '-v7.3');
fprintf('Saved wide ROI medium search to:\n%s\n', out_dir);
end

function data = generate_case_data(scene, cfg, snr_db, W_bands)
num_bands = numel(scene.freq_ratios);
z_bands = zeros(cfg.P, num_bands);
D_bands = zeros(cfg.P, numel(cfg.scan_deg), num_bands);
c_bands = zeros(cfg.P, num_bands);
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

function W = build_W_from_beams(beam_deg, M, d_over_lambda)
W = zeros(numel(beam_deg), M);
for p = 1:numel(beam_deg)
    W(p, :) = steering_vector_local(beam_deg(p), M, d_over_lambda)';
end
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

function cfg_old = make_old_cfg(cfg)
cfg_old = struct();
cfg_old.L = cfg.L;
cfg_old.d_over_lambda = cfg.d_over_lambda;
cfg_old.fit_floor = 1e-8;
cfg_old.max_alt_iters = 12;
cfg_old.max_refine_passes = 2;
cfg_old.refine_radii_deg = [1.0, 0.5, 0.25, 0.12, 0.06];
cfg_old.refine_samples = 11;
cfg_old.refine_guard_deg = 0.6;
cfg_old.angle_min_deg = cfg.scan_deg(1);
cfg_old.angle_max_deg = cfg.scan_deg(end);
cfg_old.freq_ratios = cfg.freq_ratios;
cfg_old.max_support_evals = 50000;
end
