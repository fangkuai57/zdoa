function summary = benchmark_zdoa_ongrid_k123(quick_mode)
if nargin < 1
    quick_mode = false;
end

close all;
clc;
rng(2026);

cfg.M = 20;
cfg.d_over_lambda = 0.5;
cfg.scan_deg = -18:2:18;
cfg.max_support_evals = 50000;
cfg.max_alt_iters = 12;
cfg.fit_floor = 1e-8;

k_specs = {
    struct('K', 1, 'name', 'K=1', 'source_power', 1.0, 'min_sep_deg', 0), ...
    struct('K', 2, 'name', 'K=2', 'source_power', [1.0, 1.0], 'min_sep_deg', 4), ...
    struct('K', 3, 'name', 'K=3', 'source_power', [1.0, 0.9, 0.7], 'min_sep_deg', 4)
};

if quick_mode
    num_trials = 24;
    beam_count_values = [6, 8, 10, 12, 16, 20];
    snr_values = -10:2:10;
    compression_setup.snr_db = 4;
    compression_setup.L = 40;
    snr_setup.beam_count = 20;
    snr_setup.L = 40;
else
    num_trials = 72;
    beam_count_values = [6, 8, 10, 12, 16, 20];
    snr_values = -12:2:12;
    compression_setup.snr_db = 4;
    compression_setup.L = 40;
    snr_setup.beam_count = 20;
    snr_setup.L = 40;
end

reliable_rate = 0.90;
summary.meta.algorithm = 'zdoa_ongrid_exact';
summary.meta.num_trials = num_trials;
summary.meta.reliable_rate = reliable_rate;
summary.meta.M = cfg.M;
summary.meta.scan_deg = cfg.scan_deg;
summary.meta.beam_count_values = beam_count_values;
summary.meta.compression_ratio_values = beam_count_values / cfg.M;
summary.meta.snr_values = snr_values;
summary.meta.k_values = cellfun(@(s) s.K, k_specs);

summary.compression = run_compression_benchmark(cfg, k_specs, beam_count_values, compression_setup, num_trials);
summary.snr = run_snr_benchmark(cfg, k_specs, snr_values, snr_setup, num_trials);

for idx = 1:numel(k_specs)
    k_name = matlab.lang.makeValidName(sprintf('K%d', k_specs{idx}.K));
    summary.boundary.(k_name).min_reliable_beam_count = first_reliable( ...
        summary.compression.per_k(idx).beam_count_values, ...
        summary.compression.per_k(idx).exact_rate, ...
        reliable_rate);
    summary.boundary.(k_name).min_reliable_snr_db = first_reliable( ...
        summary.snr.per_k(idx).snr_values, ...
        summary.snr.per_k(idx).exact_rate, ...
        reliable_rate);
end

fig = figure('Color', 'w', 'Position', [60, 60, 1320, 920]);
colors = lines(numel(k_specs));

subplot(2, 2, 1);
hold on;
for idx = 1:numel(k_specs)
    plot(summary.compression.per_k(idx).beam_count_values, ...
         summary.compression.per_k(idx).exact_rate, ...
         '-o', 'Color', colors(idx, :), 'LineWidth', 1.6, 'MarkerSize', 5, ...
         'DisplayName', k_specs{idx}.name);
end
grid on;
ylim([0, 1.02]);
xlabel('Front-end output size P');
ylabel('Exact recovery rate');
title(sprintf('Compression sweep (L=%d, SNR=%d dB)', compression_setup.L, compression_setup.snr_db));
legend('Location', 'southeast');

subplot(2, 2, 2);
hold on;
for idx = 1:numel(k_specs)
    plot(summary.compression.per_k(idx).beam_count_values, ...
         summary.compression.per_k(idx).mean_mse, ...
         '-o', 'Color', colors(idx, :), 'LineWidth', 1.6, 'MarkerSize', 5, ...
         'DisplayName', k_specs{idx}.name);
end
grid on;
xlabel('Front-end output size P');
ylabel('Mean squared DOA error (deg^2)');
title(sprintf('Compression sweep MSE (L=%d, SNR=%d dB)', compression_setup.L, compression_setup.snr_db));
legend('Location', 'northeast');

subplot(2, 2, 3);
hold on;
for idx = 1:numel(k_specs)
    plot(summary.snr.per_k(idx).snr_values, ...
         summary.snr.per_k(idx).exact_rate, ...
         '-o', 'Color', colors(idx, :), 'LineWidth', 1.6, 'MarkerSize', 5, ...
         'DisplayName', k_specs{idx}.name);
end
grid on;
ylim([0, 1.02]);
xlabel('SNR (dB)');
ylabel('Exact recovery rate');
title(sprintf('K sweep under fixed P=%d, L=%d', snr_setup.beam_count, snr_setup.L));
legend('Location', 'southeast');

subplot(2, 2, 4);
hold on;
for idx = 1:numel(k_specs)
    plot(summary.snr.per_k(idx).snr_values, ...
         summary.snr.per_k(idx).mean_mse, ...
         '-o', 'Color', colors(idx, :), 'LineWidth', 1.6, 'MarkerSize', 5, ...
         'DisplayName', k_specs{idx}.name);
end
grid on;
xlabel('SNR (dB)');
ylabel('Mean squared DOA error (deg^2)');
title(sprintf('K sweep MSE under fixed P=%d, L=%d', snr_setup.beam_count, snr_setup.L));
legend('Location', 'northeast');

sgtitle(sprintf('On-grid exact benchmark: compression and K=1/2/3, scan grid = %d points', numel(cfg.scan_deg)));

if quick_mode
    output_png = 'benchmark_zdoa_ongrid_k123_quick.png';
    output_mat = 'benchmark_zdoa_ongrid_k123_quick.mat';
else
    output_png = 'benchmark_zdoa_ongrid_k123.png';
    output_mat = 'benchmark_zdoa_ongrid_k123.mat';
end

exportgraphics(fig, output_png, 'Resolution', 180);
save(output_mat, 'summary', 'cfg', 'k_specs', 'compression_setup', 'snr_setup');

fprintf('\n=== On-grid K=1/2/3 benchmark summary ===\n');
fprintf('Compression sweep: L = %d, SNR = %d dB\n', compression_setup.L, compression_setup.snr_db);
fprintf('SNR sweep: fixed P = %d, L = %d\n', snr_setup.beam_count, snr_setup.L);
for idx = 1:numel(k_specs)
    k_name = sprintf('K%d', k_specs{idx}.K);
    fprintf('%s: min reliable P = %g, min reliable SNR = %g dB\n', ...
        k_name, ...
        summary.boundary.(k_name).min_reliable_beam_count, ...
        summary.boundary.(k_name).min_reliable_snr_db);
end
end

function result = run_compression_benchmark(base_cfg, k_specs, beam_count_values, setup, num_trials)
per_k = struct([]);

for idx = 1:numel(k_specs)
    exact_rate = zeros(size(beam_count_values));
    mean_mse = zeros(size(beam_count_values));

    cfg = base_cfg;
    cfg.K = k_specs{idx}.K;
    cfg.L = setup.L;
    cfg.snr_db = setup.snr_db;

    support_bank = sample_support_bank(cfg.scan_deg, cfg.K, num_trials, k_specs{idx}.min_sep_deg);
    for beam_idx = 1:numel(beam_count_values)
        beam_deg = linspace(cfg.scan_deg(1), cfg.scan_deg(end), beam_count_values(beam_idx));
        W = build_beamformer_local(cfg.M, beam_deg, cfg.d_over_lambda);
        [exact_rate(beam_idx), mean_mse(beam_idx)] = evaluate_trials(cfg, W, support_bank, k_specs{idx}.source_power);
    end

    per_k(idx).K = k_specs{idx}.K;
    per_k(idx).name = k_specs{idx}.name;
    per_k(idx).beam_count_values = beam_count_values;
    per_k(idx).compression_ratio_values = beam_count_values / base_cfg.M;
    per_k(idx).exact_rate = exact_rate;
    per_k(idx).mean_mse = mean_mse;
    per_k(idx).table = table(beam_count_values(:), (beam_count_values(:) / base_cfg.M), exact_rate(:), mean_mse(:), ...
        'VariableNames', {'beam_count', 'compression_ratio', 'exact_rate', 'mean_mse'});
end

result.per_k = per_k;
end

function result = run_snr_benchmark(base_cfg, k_specs, snr_values, setup, num_trials)
per_k = struct([]);

beam_deg = linspace(base_cfg.scan_deg(1), base_cfg.scan_deg(end), setup.beam_count);
W = build_beamformer_local(base_cfg.M, beam_deg, base_cfg.d_over_lambda);

for idx = 1:numel(k_specs)
    exact_rate = zeros(size(snr_values));
    mean_mse = zeros(size(snr_values));

    cfg = base_cfg;
    cfg.K = k_specs{idx}.K;
    cfg.L = setup.L;

    support_bank = sample_support_bank(cfg.scan_deg, cfg.K, num_trials, k_specs{idx}.min_sep_deg);
    for snr_idx = 1:numel(snr_values)
        cfg.snr_db = snr_values(snr_idx);
        [exact_rate(snr_idx), mean_mse(snr_idx)] = evaluate_trials(cfg, W, support_bank, k_specs{idx}.source_power);
    end

    per_k(idx).K = k_specs{idx}.K;
    per_k(idx).name = k_specs{idx}.name;
    per_k(idx).snr_values = snr_values;
    per_k(idx).exact_rate = exact_rate;
    per_k(idx).mean_mse = mean_mse;
    per_k(idx).table = table(snr_values(:), exact_rate(:), mean_mse(:), ...
        'VariableNames', {'snr_db', 'exact_rate', 'mean_mse'});
end

result.per_k = per_k;
result.fixed_beam_deg = beam_deg;
result.fixed_frontend_output_size = [size(W, 1), 1];
end

function [exact_rate, mean_mse] = evaluate_trials(cfg, W, support_bank, source_power)
num_trials = size(support_bank, 1);
success = false(num_trials, 1);
mse_values = zeros(num_trials, 1);

for trial = 1:num_trials
    cfg.true_doa_deg = support_bank(trial, :);
    cfg.source_power = source_power;
    X = generate_snapshots_local(cfg);
    z = mean(abs(W * X).^2, 2);
    [~, doa_est] = zdoa_ongrid_exact(z, W, cfg.scan_deg, cfg.K, cfg);
    doa_true = sort(cfg.true_doa_deg(:)).';
    doa_hat = sort(doa_est(:)).';
    success(trial) = isequal(doa_hat, doa_true);
    mean_len = min(numel(doa_hat), numel(doa_true));
    mse_values(trial) = mean((doa_hat(1:mean_len) - doa_true(1:mean_len)).^2);
end

exact_rate = mean(success);
mean_mse = mean(mse_values);
end

function support_bank = sample_support_bank(scan_deg, K, num_trials, min_sep_deg)
scan_deg = scan_deg(:).';
all_supports = nchoosek(scan_deg, K);

if K > 1 && min_sep_deg > 0
    keep = true(size(all_supports, 1), 1);
    for row = 1:size(all_supports, 1)
        if any(diff(all_supports(row, :)) < min_sep_deg)
            keep(row) = false;
        end
    end
    all_supports = all_supports(keep, :);
end

if isempty(all_supports)
    error('No valid supports found for K=%d and minimum separation %g deg.', K, min_sep_deg);
end

support_bank = zeros(num_trials, K);
choices = randi(size(all_supports, 1), num_trials, 1);
for trial = 1:num_trials
    support_bank(trial, :) = all_supports(choices(trial), :);
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
