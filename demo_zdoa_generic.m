function results = demo_zdoa_generic()
close all;
clc;
rng(42);

cfg.M = 16;
cfg.d_over_lambda = 0.5;
cfg.K = 2;
cfg.scan_deg = -40:0.25:40;
cfg.lambda_sparse = 0.04;
cfg.lambda_smooth = 0.01;
cfg.reweight_loops = 4;
cfg.reweight_epsilon = 1e-3;
cfg.max_iters = 350;
cfg.step_shrink = 0.98;
cfg.power_weight_floor = 1e-6;
cfg.refine_span_deg = 1.0;
cfg.candidate_pool_size = 10;
cfg.max_refine_candidates = 16;

beam_deg = linspace(-40, 40, 11);
W = build_beamformer_local(cfg.M, beam_deg, cfg.d_over_lambda);

cases = {
    struct('name', 'Well separated, medium SNR', 'true_doa_deg', [-22, 16], 'source_power', [1.0, 0.9], 'L', 64, 'snr_db', 2), ...
    struct('name', 'Moderately close, equal power', 'true_doa_deg', [-11, -4], 'source_power', [1.0, 1.0], 'L', 72, 'snr_db', 6), ...
    struct('name', 'Strong-weak adjacent sources', 'true_doa_deg', [-6, -1], 'source_power', [1.0, 0.35], 'L', 72, 'snr_db', 6), ...
    struct('name', 'Low SNR or limited snapshots', 'true_doa_deg', [-16, 10], 'source_power', [1.0, 0.8], 'L', 16, 'snr_db', -4)
};

results = struct([]);
fig = figure('Color', 'w', 'Position', [80, 80, 1280, 860]);

for case_idx = 1:numel(cases)
    run_cfg = cfg;
    run_cfg.true_doa_deg = cases{case_idx}.true_doa_deg;
    run_cfg.source_power = cases{case_idx}.source_power;
    run_cfg.L = cases{case_idx}.L;
    run_cfg.snr_db = cases{case_idx}.snr_db;

    X = generate_snapshots_local(run_cfg);
    Y = W * X;
    z = mean(abs(Y).^2, 2);
    [spectrum, doa_est, debug] = zdoa_generic_recover(z, W, run_cfg.scan_deg, run_cfg.K, run_cfg);

    results(case_idx).name = cases{case_idx}.name;
    results(case_idx).true_doa_deg = run_cfg.true_doa_deg;
    results(case_idx).source_power = run_cfg.source_power;
    results(case_idx).L = run_cfg.L;
    results(case_idx).snr_db = run_cfg.snr_db;
    results(case_idx).estimated_doa_deg = doa_est;
    results(case_idx).beam_output_size = size(z);
    results(case_idx).spectrum = spectrum;
    results(case_idx).debug = debug;

    fprintf('\n=== %s ===\n', cases{case_idx}.name);
    fprintf('True DOA (deg):      '); fprintf('%.2f ', run_cfg.true_doa_deg); fprintf('\n');
    fprintf('Estimated DOA (deg): '); fprintf('%.2f ', doa_est); fprintf('\n');
    fprintf('Front-end output size z: %d x 1\n', numel(z));

    subplot(2, 2, case_idx);
    plot(run_cfg.scan_deg, spectrum, 'r-', 'LineWidth', 1.4);
    hold on;
    for k = 1:numel(run_cfg.true_doa_deg)
        xline(run_cfg.true_doa_deg(k), '--k', 'LineWidth', 1.2);
    end
    for k = 1:numel(doa_est)
        xline(doa_est(k), ':b', 'LineWidth', 1.2);
    end
    grid on;
    ylim([0, 1.05]);
    xlabel('Angle (deg)');
    ylabel('Normalized spectrum');
    title(sprintf('%s\nL=%d, SNR=%d dB', cases{case_idx}.name, run_cfg.L, run_cfg.snr_db));
    legend('Recovered spectrum', 'True DOA', 'Estimated DOA', 'Location', 'northeast');
end

sgtitle(sprintf('Generic Z-domain DOA demo, front-end output size = %d x 1', size(z, 1)));
exportgraphics(fig, 'demo_zdoa_generic.png', 'Resolution', 180);
save('demo_zdoa_generic.mat', 'results', 'cfg', 'W', 'beam_deg');
end

function W = build_beamformer_local(M, beam_deg, d_over_lambda)
B = numel(beam_deg);
W = zeros(B, M);
for b = 1:B
    W(b, :) = steering_vector_local(beam_deg(b), M, d_over_lambda)';
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

function a = steering_vector_local(theta_deg, M, d_over_lambda)
n = (0:M-1).';
a = exp(1i * 2 * pi * d_over_lambda * n * sind(theta_deg));
a = a / sqrt(M);
end
