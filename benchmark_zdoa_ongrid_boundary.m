function summary = benchmark_zdoa_ongrid_boundary(quick_mode)
if nargin < 1
    quick_mode = false;
end

close all;
clc;
rng(2026);

cfg.M = 16;
cfg.d_over_lambda = 0.5;
cfg.K = 2;
cfg.scan_deg = -24:1:24;
cfg.max_support_evals = 50000;
cfg.max_alt_iters = 12;
cfg.fit_floor = 1e-8;
beam_deg = -24:4:24;
W = build_beamformer_local(cfg.M, beam_deg, cfg.d_over_lambda);

if quick_mode
    num_trials = 16;
    snr_values = -4:2:14;
    snapshot_values = [8, 12, 16, 24, 32, 48, 64];
    separation_values = [2, 3, 4, 5, 6, 8, 10, 12];
    weak_ratio_values = [1.0, 0.8, 0.6, 0.4, 0.3, 0.2];
else
    num_trials = 60;
    snr_values = -6:2:16;
    snapshot_values = [6, 8, 12, 16, 24, 32, 48, 64, 96];
    separation_values = [2, 3, 4, 5, 6, 8, 10, 12, 14];
    weak_ratio_values = [1.0, 0.8, 0.6, 0.4, 0.3, 0.2, 0.1];
end

reliable_rate = 0.90;

eval_cfg.base_center_deg = 0;

summary.meta.num_trials = num_trials;
summary.meta.reliable_rate = reliable_rate;
summary.meta.frontend_output_size = [size(W, 1), 1];
summary.meta.scan_deg = cfg.scan_deg;
summary.meta.beam_deg = beam_deg;
summary.meta.algorithm = 'zdoa_ongrid_exact';

snr_setup.L = 64;
snr_setup.sep_deg = 8;
snr_setup.weak_ratio = 1.0;
summary.snr = run_sweep(cfg, W, snr_values, num_trials, 'snr', snr_setup);

snapshot_setup.snr_db = 8;
snapshot_setup.sep_deg = 8;
snapshot_setup.weak_ratio = 1.0;
summary.snapshots = run_sweep(cfg, W, snapshot_values, num_trials, 'snapshots', snapshot_setup);

separation_setup.snr_db = 10;
separation_setup.L = 96;
separation_setup.weak_ratio = 1.0;
summary.separation = run_sweep(cfg, W, separation_values, num_trials, 'separation', separation_setup);

weak_setup.snr_db = 10;
weak_setup.L = 96;
weak_setup.sep_deg = 8;
summary.weak_ratio = run_sweep(cfg, W, weak_ratio_values, num_trials, 'weak_ratio', weak_setup);

summary.boundary.min_reliable_snr_db = first_reliable(summary.snr.x_values, summary.snr.exact_rate, reliable_rate);
summary.boundary.min_reliable_snapshots = first_reliable(summary.snapshots.x_values, summary.snapshots.exact_rate, reliable_rate);
summary.boundary.min_reliable_separation_deg = first_reliable(summary.separation.x_values, summary.separation.exact_rate, reliable_rate);
summary.boundary.min_reliable_weak_ratio = last_reliable(summary.weak_ratio.x_values, summary.weak_ratio.exact_rate, reliable_rate);

fig = figure('Color', 'w', 'Position', [80, 80, 1200, 860]);
subplot(2, 2, 1);
plot(summary.snr.x_values, summary.snr.exact_rate, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
ylim([0, 1.02]);
xlabel('SNR (dB)');
ylabel('Exact recovery rate');
title(sprintf('SNR sweep (L=%d, sep=%d^o)', snr_setup.L, snr_setup.sep_deg));

subplot(2, 2, 2);
plot(summary.snapshots.x_values, summary.snapshots.exact_rate, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
ylim([0, 1.02]);
xlabel('Number of snapshots L');
ylabel('Exact recovery rate');
title(sprintf('Snapshot sweep (SNR=%d dB, sep=%d^o)', snapshot_setup.snr_db, snapshot_setup.sep_deg));

subplot(2, 2, 3);
plot(summary.separation.x_values, summary.separation.exact_rate, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
ylim([0, 1.02]);
xlabel('Separation (deg)');
ylabel('Exact recovery rate');
title(sprintf('Separation sweep (SNR=%d dB, L=%d)', separation_setup.snr_db, separation_setup.L));

subplot(2, 2, 4);
plot(summary.weak_ratio.x_values, summary.weak_ratio.exact_rate, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on;
ylim([0, 1.02]);
xlabel('Weak/strong power ratio');
ylabel('Exact recovery rate');
title(sprintf('Power imbalance sweep (SNR=%d dB, L=%d)', weak_setup.snr_db, weak_setup.L));

sgtitle(sprintf('On-grid boundary benchmark, front-end output size = %d x 1', size(W, 1)));

if quick_mode
    output_png = 'benchmark_zdoa_ongrid_boundary_quick.png';
    output_mat = 'benchmark_zdoa_ongrid_boundary_quick.mat';
else
    output_png = 'benchmark_zdoa_ongrid_boundary.png';
    output_mat = 'benchmark_zdoa_ongrid_boundary.mat';
end

exportgraphics(fig, output_png, 'Resolution', 180);
save(output_mat, 'summary', 'cfg', 'W', 'beam_deg');

fprintf('\n=== On-grid boundary summary ===\n');
fprintf('Front-end output size: %d x 1\n', size(W, 1));
fprintf('Min reliable SNR (>= %.2f exact rate): %g dB\n', reliable_rate, summary.boundary.min_reliable_snr_db);
fprintf('Min reliable snapshots (>= %.2f exact rate): %g\n', reliable_rate, summary.boundary.min_reliable_snapshots);
fprintf('Min reliable separation (>= %.2f exact rate): %g deg\n', reliable_rate, summary.boundary.min_reliable_separation_deg);
fprintf('Smallest reliable weak ratio (>= %.2f exact rate): %g\n', reliable_rate, summary.boundary.min_reliable_weak_ratio);
end

function result = run_sweep(base_cfg, W, x_values, num_trials, mode_name, setup)
exact_rate = zeros(size(x_values));
mean_mse = zeros(size(x_values));

for idx = 1:numel(x_values)
    cfg = base_cfg;
    switch mode_name
        case 'snr'
            cfg.snr_db = x_values(idx);
            cfg.L = setup.L;
            sep_deg = setup.sep_deg;
            weak_ratio = setup.weak_ratio;
        case 'snapshots'
            cfg.snr_db = setup.snr_db;
            cfg.L = x_values(idx);
            sep_deg = setup.sep_deg;
            weak_ratio = setup.weak_ratio;
        case 'separation'
            cfg.snr_db = setup.snr_db;
            cfg.L = setup.L;
            sep_deg = x_values(idx);
            weak_ratio = setup.weak_ratio;
        case 'weak_ratio'
            cfg.snr_db = setup.snr_db;
            cfg.L = setup.L;
            sep_deg = setup.sep_deg;
            weak_ratio = x_values(idx);
        otherwise
            error('Unknown sweep mode: %s', mode_name);
    end

    success = false(num_trials, 1);
    mse_values = zeros(num_trials, 1);
    pair_bank = sample_pair_bank(cfg.scan_deg, num_trials, sep_deg);

    for trial = 1:num_trials
        cfg.true_doa_deg = pair_bank(trial, :);
        cfg.source_power = [1.0, weak_ratio];
        X = generate_snapshots_local(cfg);
        z = mean(abs(W * X).^2, 2);
        [~, doa_est] = zdoa_ongrid_exact(z, W, cfg.scan_deg, cfg.K, cfg);
        doa_true = sort(cfg.true_doa_deg(:)).';
        doa_hat = sort(doa_est(:)).';
        success(trial) = isequal(doa_hat, doa_true);
        mse_values(trial) = mean((doa_hat - doa_true).^2);
    end

    exact_rate(idx) = mean(success);
    mean_mse(idx) = mean(mse_values);
end

result.x_values = x_values;
result.exact_rate = exact_rate;
result.mean_mse = mean_mse;
result.table = table(x_values(:), exact_rate(:), mean_mse(:), 'VariableNames', {'x_value', 'exact_rate', 'mean_mse'});
end

function pair_bank = sample_pair_bank(scan_deg, num_trials, sep_deg)
scan_deg = scan_deg(:).';
pair_bank = zeros(num_trials, 2);
valid_pairs = [];
for idx1 = 1:numel(scan_deg)
    for idx2 = idx1 + 1:numel(scan_deg)
        if abs(scan_deg(idx2) - scan_deg(idx1)) == sep_deg
            valid_pairs = [valid_pairs; scan_deg(idx1), scan_deg(idx2)]; %#ok<AGROW>
        end
    end
end
if isempty(valid_pairs)
    error('No valid on-grid pairs found for separation %g deg.', sep_deg);
end
for trial = 1:num_trials
    choice = randi(size(valid_pairs, 1));
    pair_bank(trial, :) = valid_pairs(choice, :);
end
end

function X = generate_snapshots_local(cfg)
M = cfg.M;
L = cfg.L;
K = cfg.K;
A = zeros(M, K);
for k = 1:K
    A(:, k) = steering_vector_local(cfg.true_doa_deg(k), M, cfg.d_over_lambda);
end
S = (randn(K, L) + 1i * randn(K, L)) / sqrt(2);
for k = 1:K
    S(k, :) = sqrt(cfg.source_power(k)) * S(k, :);
end
signal_var = mean(abs(A * S).^2, 'all');
noise_var = signal_var / (10^(cfg.snr_db / 10));
N = sqrt(noise_var / 2) * (randn(M, L) + 1i * randn(M, L));
X = A * S + N;
end

function W = build_beamformer_local(M, beam_deg, d_over_lambda)
B = numel(beam_deg);
W = zeros(B, M);
for b = 1:B
    W(b, :) = steering_vector_local(beam_deg(b), M, d_over_lambda)';
end
end

function a = steering_vector_local(theta_deg, M, d_over_lambda)
n = (0:M-1).';
a = exp(1i * 2 * pi * d_over_lambda * n * sind(theta_deg));
a = a / sqrt(M);
end

function value = first_reliable(x_values, rates, threshold)
idx = find(rates >= threshold, 1, 'first');
if isempty(idx)
    value = NaN;
else
    value = x_values(idx);
end
end

function value = last_reliable(x_values, rates, threshold)
idx = find(rates >= threshold, 1, 'last');
if isempty(idx)
    value = NaN;
else
    value = x_values(idx);
end
end
