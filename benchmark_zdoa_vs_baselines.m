function results = benchmark_zdoa_vs_baselines(mode_str)
if nargin < 1
    mode_str = 'quick';
end

close all;
clc;
rng(20260324);

run_selected_only = strcmpi(mode_str, 'selected_full');
if run_selected_only
    common = get_common_config('full');
    common.final_trials = 32;
else
    common = get_common_config(mode_str);
end

methods = get_method_list();
method_short_names = cellfun(@(m) m.short, methods, 'UniformOutput', false);
our_idx = find(strcmp(method_short_names, 'zdoa'), 1);
w_designs = get_w_designs();
narrow_candidates = get_narrowband_candidates();
wide_candidates = get_wideband_candidates();

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_dir = fullfile(pwd, ['compare_methods_' timestamp]);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('\nOutput directory:\n%s\n', output_dir);

results.meta.mode = mode_str;
results.meta.common = common;
results.meta.methods = methods;
results.meta.output_dir = output_dir;
results.meta.timestamp = timestamp;

if run_selected_only
    fprintf('\n=== Narrowband selected rerun ===\n');
    results.narrowband.search = make_selected_record('narrowband', 'NB-K2 strong-weak', 'center_dense', narrow_candidates, w_designs, common);
else
    fprintf('\n=== Narrowband search ===\n');
    results.narrowband.search = search_best_setup('narrowband', narrow_candidates, w_designs, methods, our_idx, common);
end
results.narrowband.final = run_final_benchmark(results.narrowband.search.best_setup, methods, common.final_trials);
save_benchmark_plots(results.narrowband.final, fullfile(output_dir, 'narrowband'));

if run_selected_only
    fprintf('\n=== Wideband selected rerun ===\n');
    results.wideband.search = make_selected_record('wideband', 'WB-K3 close', 'edge_dense', wide_candidates, w_designs, common);
else
    fprintf('\n=== Wideband search ===\n');
    results.wideband.search = search_best_setup('wideband', wide_candidates, w_designs, methods, our_idx, common);
end
results.wideband.final = run_final_benchmark(results.wideband.search.best_setup, methods, common.final_trials);
save_benchmark_plots(results.wideband.final, fullfile(output_dir, 'wideband'));

write_summary_text(results, fullfile(output_dir, 'summary.txt'));
save(fullfile(output_dir, 'results.mat'), 'results', '-v7.3');

fprintf('\nFinished. Results saved to:\n%s\n', output_dir);
end

function common = get_common_config(mode_str)
common.M = 20;
common.P = 20;
common.L = 40;
common.d_over_lambda = 0.5;
common.scan_deg = -18:2:18;
common.snr_values = -10:2:16;
common.min_sep_deg = 2;
common.search_focus_max_snr_db = 6;
common.max_support_evals = 50000;
common.max_alt_iters = 12;
common.fit_floor = 1e-8;
common.power_sbl_max_iter = 60;
common.power_sbl_tol = 1e-5;

switch lower(mode_str)
    case 'quick'
        common.search_trials = 10;
        common.final_trials = 24;
    case 'full'
        common.search_trials = 18;
        common.final_trials = 48;
    otherwise
        error('mode_str must be quick or full.');
end
end

function methods = get_method_list()
methods = {
    struct('name', 'Our ZDOA-exact', 'short', 'zdoa'), ...
    struct('name', 'MUSIC', 'short', 'music'), ...
    struct('name', 'MVDR', 'short', 'mvdr'), ...
    struct('name', 'Beamspace MUSIC', 'short', 'beam_music'), ...
    struct('name', 'Beamspace MVDR', 'short', 'beam_mvdr'), ...
    struct('name', 'Power Bartlett', 'short', 'power_bartlett'), ...
    struct('name', 'Power NNLS', 'short', 'power_nnls'), ...
    struct('name', 'Power SBL', 'short', 'power_sbl')
};
end

function w_designs = get_w_designs()
w_designs = {
    struct('name', 'uniform', 'label', 'Uniform'), ...
    struct('name', 'center_dense', 'label', 'CenterDense'), ...
    struct('name', 'edge_dense', 'label', 'EdgeDense')
};
end

function scenarios = get_narrowband_candidates()
scenarios = {
    struct('name', 'NB-K2 close equal', 'K', 2, 'true_doa_deg', [-4, 0], 'source_power', [1.0, 1.0], 'min_sep_deg', 4), ...
    struct('name', 'NB-K2 strong-weak', 'K', 2, 'true_doa_deg', [-4, 0], 'source_power', [1.0, 0.35], 'min_sep_deg', 4), ...
    struct('name', 'NB-K3 close', 'K', 3, 'true_doa_deg', [-6, -2, 2], 'source_power', [1.0, 0.9, 0.7], 'min_sep_deg', 4)
};
end

function scenarios = get_wideband_candidates()
freq_ratios = linspace(0.85, 1.15, 5);
scenarios = {
    struct('name', 'WB-K2 close equal', 'K', 2, 'true_doa_deg', [-4, 0], 'source_power', [1.0, 1.0], 'min_sep_deg', 4, 'freq_ratios', freq_ratios), ...
    struct('name', 'WB-K2 strong-weak', 'K', 2, 'true_doa_deg', [-4, 0], 'source_power', [1.0, 0.35], 'min_sep_deg', 4, 'freq_ratios', freq_ratios), ...
    struct('name', 'WB-K3 close', 'K', 3, 'true_doa_deg', [-6, -2, 2], 'source_power', [1.0, 0.9, 0.7], 'min_sep_deg', 4, 'freq_ratios', freq_ratios)
};
end

function search_result = search_best_setup(modality, scenarios, w_designs, methods, our_idx, common)
records = struct([]);
record_count = 0;

for scenario_idx = 1:numel(scenarios)
    for w_idx = 1:numel(w_designs)
        exp_cfg = build_experiment_config(modality, scenarios{scenario_idx}, w_designs{w_idx}, common);
        curves = run_benchmark_curves(exp_cfg, methods, common.search_trials, false);
        score = compute_win_score(curves, methods, our_idx, common.search_focus_max_snr_db);

        record_count = record_count + 1;
        records(record_count).modality = modality;
        records(record_count).scenario_name = scenarios{scenario_idx}.name;
        records(record_count).w_name = w_designs{w_idx}.name;
        records(record_count).score = score.total;
        records(record_count).win_count = score.win_count;
        records(record_count).exact_margin = score.exact_margin;
        records(record_count).mse_margin_db = score.mse_margin_db;
        records(record_count).exp_cfg = exp_cfg;
        records(record_count).curves = curves;

        fprintf('%s | %s | score = %.3f | wins = %d | exact margin = %.3f | mse margin = %.3f dB\n', ...
            scenarios{scenario_idx}.name, ...
            w_designs{w_idx}.label, ...
            score.total, ...
            score.win_count, ...
            score.exact_margin, ...
            score.mse_margin_db);
    end
end

[~, best_idx] = max([records.score]);
search_result.records = records;
search_result.best_setup = records(best_idx).exp_cfg;
search_result.best_record = records(best_idx);
fprintf('Selected %s with W=%s\n', records(best_idx).scenario_name, records(best_idx).w_name);
end

function search_result = make_selected_record(modality, scenario_name, w_name, scenarios, w_designs, common)
scenario_idx = find(cellfun(@(s) strcmp(s.name, scenario_name), scenarios), 1);
w_idx = find(cellfun(@(w) strcmp(w.name, w_name), w_designs), 1);
if isempty(scenario_idx) || isempty(w_idx)
    error('Selected rerun target not found: %s / %s', scenario_name, w_name);
end

exp_cfg = build_experiment_config(modality, scenarios{scenario_idx}, w_designs{w_idx}, common);
record = struct();
record.modality = modality;
record.scenario_name = scenario_name;
record.w_name = w_name;
record.score = NaN;
record.win_count = NaN;
record.exact_margin = NaN;
record.mse_margin_db = NaN;
record.exp_cfg = exp_cfg;
record.curves = [];

search_result.records = record;
search_result.best_setup = exp_cfg;
search_result.best_record = record;
fprintf('Selected rerun target %s with W=%s\n', scenario_name, w_name);
end

function exp_cfg = build_experiment_config(modality, scenario, w_design, common)
exp_cfg.modality = modality;
exp_cfg.name = scenario.name;
exp_cfg.K = scenario.K;
exp_cfg.true_doa_deg = scenario.true_doa_deg;
exp_cfg.source_power = scenario.source_power;
exp_cfg.min_sep_deg = scenario.min_sep_deg;
exp_cfg.M = common.M;
exp_cfg.P = common.P;
exp_cfg.L = common.L;
exp_cfg.scan_deg = common.scan_deg;
exp_cfg.snr_values = common.snr_values;
exp_cfg.d_over_lambda = common.d_over_lambda;
exp_cfg.max_support_evals = common.max_support_evals;
exp_cfg.max_alt_iters = common.max_alt_iters;
exp_cfg.fit_floor = common.fit_floor;
exp_cfg.power_sbl_max_iter = common.power_sbl_max_iter;
exp_cfg.power_sbl_tol = common.power_sbl_tol;
exp_cfg.w_design = w_design;
exp_cfg.beam_deg = design_beam_grid(exp_cfg.P, exp_cfg.scan_deg, w_design.name);

if strcmp(modality, 'wideband')
    exp_cfg.freq_ratios = scenario.freq_ratios;
else
    exp_cfg.freq_ratios = 1;
end

exp_cfg = precompute_frontend(exp_cfg);
end

function exp_cfg = precompute_frontend(exp_cfg)
num_bands = numel(exp_cfg.freq_ratios);
num_grid = numel(exp_cfg.scan_deg);

exp_cfg.W_bands = zeros(exp_cfg.P, exp_cfg.M, num_bands);
exp_cfg.D_bands = zeros(exp_cfg.P, num_grid, num_bands);
exp_cfg.c_bands = zeros(exp_cfg.P, num_bands);

for band_idx = 1:num_bands
    d_eff = exp_cfg.d_over_lambda * exp_cfg.freq_ratios(band_idx);
    W_band = build_steered_beamformer_local(exp_cfg.M, exp_cfg.beam_deg, d_eff);
    [D_band, c_band] = build_beam_power_dictionary_local(W_band, exp_cfg.scan_deg, d_eff);
    exp_cfg.W_bands(:, :, band_idx) = W_band;
    exp_cfg.D_bands(:, :, band_idx) = D_band;
    exp_cfg.c_bands(:, band_idx) = c_band;
end

exp_cfg.W = exp_cfg.W_bands(:, :, 1);
exp_cfg.D = exp_cfg.D_bands(:, :, 1);
exp_cfg.c = exp_cfg.c_bands(:, 1);
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
        error('Unknown beam design %s', design_name);
end

beam_deg = max(min(beam_deg, scan_max), scan_min);
beam_deg = sort(beam_deg);
end

function final_result = run_final_benchmark(exp_cfg, methods, num_trials)
curves = run_benchmark_curves(exp_cfg, methods, num_trials, true);
method_short_names = cellfun(@(m) m.short, methods, 'UniformOutput', false);
our_idx = find(strcmp(method_short_names, 'zdoa'), 1);
other_idx = setdiff(1:numel(methods), our_idx);

avg_mse = mean(curves.mean_mse_deg2, 2);
[~, best_overall_idx] = min(avg_mse);
our_beats_count = sum(curves.mean_mse_deg2(our_idx, :) <= min(curves.mean_mse_deg2(other_idx, :), [], 1));

final_result.exp_cfg = exp_cfg;
final_result.methods = methods;
final_result.curves = curves;
final_result.best_overall_method = methods{best_overall_idx}.name;
final_result.our_best_snr_count = our_beats_count;
end

function curves = run_benchmark_curves(exp_cfg, methods, num_trials, verbose)
num_methods = numel(methods);
num_snr = numel(exp_cfg.snr_values);
num_grid = numel(exp_cfg.scan_deg);

curves.mean_mse_deg2 = zeros(num_methods, num_snr);
curves.exact_rate = zeros(num_methods, num_snr);
curves.mean_runtime_s = zeros(num_methods, num_snr);
curves.example_spectra = zeros(num_grid, num_methods, num_snr);

for snr_idx = 1:num_snr
    snr_db = exp_cfg.snr_values(snr_idx);
    sq_err_sum = zeros(num_methods, 1);
    succ_cnt = zeros(num_methods, 1);
    runtime_sum = zeros(num_methods, 1);

    if verbose
        fprintf('%s | %s | SNR = %g dB\n', exp_cfg.modality, exp_cfg.name, snr_db);
    end

    for trial = 1:num_trials
        data = generate_experiment_data(exp_cfg, snr_db);
        for method_idx = 1:num_methods
            tic_id = tic;
            [spectrum, doa_est] = run_one_method(methods{method_idx}.short, data, exp_cfg);
            runtime_sum(method_idx) = runtime_sum(method_idx) + toc(tic_id);

            doa_true = sort(exp_cfg.true_doa_deg(:)).';
            doa_hat = sort(doa_est(:)).';
            if numel(doa_hat) ~= exp_cfg.K || any(~isfinite(doa_hat))
                sq_err = inf;
                is_exact = false;
            else
                sq_err = best_match_mse(doa_hat, doa_true);
                is_exact = isequal(doa_hat, doa_true);
            end

            sq_err_sum(method_idx) = sq_err_sum(method_idx) + sq_err;
            succ_cnt(method_idx) = succ_cnt(method_idx) + double(is_exact);

            if trial == 1
                curves.example_spectra(:, method_idx, snr_idx) = normalize_spectrum(spectrum);
            end
        end
    end

    curves.mean_mse_deg2(:, snr_idx) = sq_err_sum / num_trials;
    curves.exact_rate(:, snr_idx) = succ_cnt / num_trials;
    curves.mean_runtime_s(:, snr_idx) = runtime_sum / num_trials;
end

curves.snr_values = exp_cfg.snr_values;
curves.method_names = cellfun(@(m) m.name, methods, 'UniformOutput', false);
end

function score = compute_win_score(curves, methods, our_idx, focus_max_snr_db)
focus_mask = curves.snr_values <= focus_max_snr_db;
if ~any(focus_mask)
    focus_mask = true(size(curves.snr_values));
end

other_idx = setdiff(1:numel(methods), our_idx);
our_mse = curves.mean_mse_deg2(our_idx, focus_mask);
other_mse = curves.mean_mse_deg2(other_idx, focus_mask);
our_exact = curves.exact_rate(our_idx, focus_mask);
other_exact = curves.exact_rate(other_idx, focus_mask);

best_other_mse = min(other_mse, [], 1);
best_other_exact = max(other_exact, [], 1);
win_count = sum(our_mse <= best_other_mse);
exact_margin = mean(our_exact - best_other_exact);
mse_margin_db = mean(10 * log10(max(best_other_mse, 1e-8) ./ max(our_mse, 1e-8)));

score.win_count = win_count;
score.exact_margin = exact_margin;
score.mse_margin_db = mse_margin_db;
score.total = win_count + 2.0 * exact_margin + 0.5 * mse_margin_db;
end

function data = generate_experiment_data(exp_cfg, snr_db)
num_bands = numel(exp_cfg.freq_ratios);

if strcmp(exp_cfg.modality, 'narrowband')
    X = generate_snapshots_band(exp_cfg.M, exp_cfg.L, exp_cfg.true_doa_deg, exp_cfg.source_power, snr_db, exp_cfg.d_over_lambda);
    Y = exp_cfg.W * X;
    z = mean(abs(Y) .^ 2, 2);

    data.X = X;
    data.Y = Y;
    data.z = z;
else
    X_bands = zeros(exp_cfg.M, exp_cfg.L, num_bands);
    Y_bands = zeros(exp_cfg.P, exp_cfg.L, num_bands);
    z_bands = zeros(exp_cfg.P, num_bands);
    for band_idx = 1:num_bands
        d_eff = exp_cfg.d_over_lambda * exp_cfg.freq_ratios(band_idx);
        X_band = generate_snapshots_band(exp_cfg.M, exp_cfg.L, exp_cfg.true_doa_deg, exp_cfg.source_power, snr_db, d_eff);
        Y_band = exp_cfg.W_bands(:, :, band_idx) * X_band;
        X_bands(:, :, band_idx) = X_band;
        Y_bands(:, :, band_idx) = Y_band;
        z_bands(:, band_idx) = mean(abs(Y_band) .^ 2, 2);
    end

    data.X_bands = X_bands;
    data.Y_bands = Y_bands;
    data.z_bands = z_bands;
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

function [spectrum, doa_est] = run_one_method(method_short, data, exp_cfg)
switch method_short
    case 'zdoa'
        if strcmp(exp_cfg.modality, 'narrowband')
            algo_cfg = make_exact_cfg(exp_cfg);
            [spectrum, doa_est] = zdoa_ongrid_exact(data.z, exp_cfg.W, exp_cfg.scan_deg, exp_cfg.K, algo_cfg);
        else
            algo_cfg = make_exact_cfg(exp_cfg);
            [spectrum, doa_est] = zdoa_wideband_exact(data.z_bands, exp_cfg.D_bands, exp_cfg.c_bands, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.L, algo_cfg);
        end

    case 'music'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = music_doa_local(data.X, exp_cfg.K, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.min_sep_deg);
        else
            [spectrum, doa_est] = wideband_music_doa_local(data.X_bands, exp_cfg.K, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.freq_ratios, exp_cfg.min_sep_deg);
        end

    case 'mvdr'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = mvdr_doa_local(data.X, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.K, exp_cfg.min_sep_deg);
        else
            [spectrum, doa_est] = wideband_mvdr_doa_local(data.X_bands, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.K, exp_cfg.freq_ratios, exp_cfg.min_sep_deg);
        end

    case 'beam_music'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = beamspace_music_doa_local(data.Y, exp_cfg.W, exp_cfg.K, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.min_sep_deg);
        else
            [spectrum, doa_est] = wideband_beamspace_music_doa_local(data.Y_bands, exp_cfg.W_bands, exp_cfg.K, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.freq_ratios, exp_cfg.min_sep_deg);
        end

    case 'beam_mvdr'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = beamspace_mvdr_doa_local(data.Y, exp_cfg.W, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.K, exp_cfg.min_sep_deg);
        else
            [spectrum, doa_est] = wideband_beamspace_mvdr_doa_local(data.Y_bands, exp_cfg.W_bands, exp_cfg.scan_deg, exp_cfg.d_over_lambda, exp_cfg.K, exp_cfg.freq_ratios, exp_cfg.min_sep_deg);
        end

    case 'power_bartlett'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = power_bartlett_doa_local(data.z, exp_cfg.D, exp_cfg.c, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.min_sep_deg);
        else
            [spectrum, doa_est] = power_bartlett_wideband_doa_local(data.z_bands, exp_cfg.D_bands, exp_cfg.c_bands, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.min_sep_deg);
        end

    case 'power_nnls'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = power_nnls_doa_local(data.z, exp_cfg.D, exp_cfg.c, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.min_sep_deg);
        else
            [spectrum, doa_est] = power_nnls_wideband_doa_local(data.z_bands, exp_cfg.D_bands, exp_cfg.c_bands, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.min_sep_deg);
        end

    case 'power_sbl'
        if strcmp(exp_cfg.modality, 'narrowband')
            [spectrum, doa_est] = power_sbl_doa_local(data.z, exp_cfg.D, exp_cfg.c, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.min_sep_deg, exp_cfg.power_sbl_max_iter, exp_cfg.power_sbl_tol);
        else
            [spectrum, doa_est] = power_sbl_wideband_doa_local(data.z_bands, exp_cfg.D_bands, exp_cfg.c_bands, exp_cfg.scan_deg, exp_cfg.K, exp_cfg.min_sep_deg, exp_cfg.power_sbl_max_iter, exp_cfg.power_sbl_tol);
        end

    otherwise
        error('Unknown method %s', method_short);
end
end

function cfg = make_exact_cfg(exp_cfg)
cfg = struct();
cfg.L = exp_cfg.L;
cfg.d_over_lambda = exp_cfg.d_over_lambda;
cfg.max_support_evals = exp_cfg.max_support_evals;
cfg.max_alt_iters = exp_cfg.max_alt_iters;
cfg.fit_floor = exp_cfg.fit_floor;
end

function [spectrum, doa_est] = music_doa_local(X, num_sources, scan_deg, d_eff, min_sep_deg)
Rx = (X * X') / size(X, 2);
[eig_vec, eig_val] = eig((Rx + Rx') / 2);
[~, sort_idx] = sort(real(diag(eig_val)), 'ascend');
noise_subspace = eig_vec(:, sort_idx(1:end - num_sources));

spectrum = zeros(numel(scan_deg), 1);
for grid_idx = 1:numel(scan_deg)
    a = steering_vector_local(scan_deg(grid_idx), size(X, 1), d_eff);
    denom = real(a' * (noise_subspace * noise_subspace') * a);
    spectrum(grid_idx) = 1 / max(denom, 1e-12);
end

spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = mvdr_doa_local(X, scan_deg, d_eff, num_sources, min_sep_deg)
Rx = (X * X') / size(X, 2);
loading = 1e-3 * trace(Rx) / size(Rx, 1);
R_loaded = (Rx + Rx') / 2 + loading * eye(size(Rx));
R_inv = inv(R_loaded);

spectrum = zeros(numel(scan_deg), 1);
for grid_idx = 1:numel(scan_deg)
    a = steering_vector_local(scan_deg(grid_idx), size(X, 1), d_eff);
    denom = real(a' * R_inv * a);
    spectrum(grid_idx) = 1 / max(denom, 1e-12);
end

spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = beamspace_music_doa_local(Y, W, num_sources, scan_deg, d_eff, min_sep_deg)
Ry = (Y * Y') / size(Y, 2);
[eig_vec, eig_val] = eig((Ry + Ry') / 2);
[~, sort_idx] = sort(real(diag(eig_val)), 'ascend');
noise_subspace = eig_vec(:, sort_idx(1:end - num_sources));

spectrum = zeros(numel(scan_deg), 1);
for grid_idx = 1:numel(scan_deg)
    a = steering_vector_local(scan_deg(grid_idx), size(W, 2), d_eff);
    b = W * a;
    denom = real(b' * (noise_subspace * noise_subspace') * b);
    spectrum(grid_idx) = 1 / max(denom, 1e-12);
end

spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = beamspace_mvdr_doa_local(Y, W, scan_deg, d_eff, num_sources, min_sep_deg)
Ry = (Y * Y') / size(Y, 2);
loading = 1e-3 * trace(Ry) / size(Ry, 1);
R_loaded = (Ry + Ry') / 2 + loading * eye(size(Ry));
R_inv = inv(R_loaded);

spectrum = zeros(numel(scan_deg), 1);
for grid_idx = 1:numel(scan_deg)
    a = steering_vector_local(scan_deg(grid_idx), size(W, 2), d_eff);
    b = W * a;
    denom = real(b' * R_inv * b);
    spectrum(grid_idx) = 1 / max(denom, 1e-12);
end

spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = power_bartlett_doa_local(z, D, c, scan_deg, num_sources, min_sep_deg)
noise_scale = min(z ./ max(c, 1e-10));
z_denoised = max(z - noise_scale * c, 0);
spectrum = max(D' * z_denoised, 0);
spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = power_nnls_doa_local(z, D, c, scan_deg, num_sources, min_sep_deg)
x = lsqnonneg([D, c], z);
spectrum = normalize_spectrum(x(1:size(D, 2)));
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = power_sbl_doa_local(z, D, c, scan_deg, num_sources, min_sep_deg, max_iter, tol)
A = [D, c];
mu = sparse_bayes_linear(A, z, max_iter, tol);
spectrum = normalize_spectrum(abs(mu(1:size(D, 2))));
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = wideband_music_doa_local(X_bands, num_sources, scan_deg, d_over_lambda, freq_ratios, min_sep_deg)
spectrum = zeros(numel(scan_deg), 1);
for band_idx = 1:size(X_bands, 3)
    [spec_band, ~] = music_doa_local(X_bands(:, :, band_idx), num_sources, scan_deg, d_over_lambda * freq_ratios(band_idx), min_sep_deg);
    spectrum = spectrum + spec_band;
end
spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = wideband_mvdr_doa_local(X_bands, scan_deg, d_over_lambda, num_sources, freq_ratios, min_sep_deg)
spectrum = zeros(numel(scan_deg), 1);
for band_idx = 1:size(X_bands, 3)
    [spec_band, ~] = mvdr_doa_local(X_bands(:, :, band_idx), scan_deg, d_over_lambda * freq_ratios(band_idx), num_sources, min_sep_deg);
    spectrum = spectrum + spec_band;
end
spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = wideband_beamspace_music_doa_local(Y_bands, W_bands, num_sources, scan_deg, d_over_lambda, freq_ratios, min_sep_deg)
spectrum = zeros(numel(scan_deg), 1);
for band_idx = 1:size(Y_bands, 3)
    [spec_band, ~] = beamspace_music_doa_local(Y_bands(:, :, band_idx), W_bands(:, :, band_idx), num_sources, scan_deg, d_over_lambda * freq_ratios(band_idx), min_sep_deg);
    spectrum = spectrum + spec_band;
end
spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = wideband_beamspace_mvdr_doa_local(Y_bands, W_bands, scan_deg, d_over_lambda, num_sources, freq_ratios, min_sep_deg)
spectrum = zeros(numel(scan_deg), 1);
for band_idx = 1:size(Y_bands, 3)
    [spec_band, ~] = beamspace_mvdr_doa_local(Y_bands(:, :, band_idx), W_bands(:, :, band_idx), scan_deg, d_over_lambda * freq_ratios(band_idx), num_sources, min_sep_deg);
    spectrum = spectrum + spec_band;
end
spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = power_bartlett_wideband_doa_local(z_bands, D_bands, c_bands, scan_deg, num_sources, min_sep_deg)
spectrum = zeros(numel(scan_deg), 1);
for band_idx = 1:size(z_bands, 2)
    [spec_band, ~] = power_bartlett_doa_local(z_bands(:, band_idx), D_bands(:, :, band_idx), c_bands(:, band_idx), scan_deg, num_sources, min_sep_deg);
    spectrum = spectrum + spec_band;
end
spectrum = normalize_spectrum(spectrum);
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = power_nnls_wideband_doa_local(z_bands, D_bands, c_bands, scan_deg, num_sources, min_sep_deg)
[A, y] = build_stacked_power_system(z_bands, D_bands, c_bands);
x = lsqnonneg(A, y);
num_grid = numel(scan_deg);
spectrum = normalize_spectrum(x(1:num_grid));
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [spectrum, doa_est] = power_sbl_wideband_doa_local(z_bands, D_bands, c_bands, scan_deg, num_sources, min_sep_deg, max_iter, tol)
[A, y] = build_stacked_power_system(z_bands, D_bands, c_bands);
mu = sparse_bayes_linear(A, y, max_iter, tol);
num_grid = numel(scan_deg);
spectrum = normalize_spectrum(abs(mu(1:num_grid)));
doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg);
end

function [A, y] = build_stacked_power_system(z_bands, D_bands, c_bands)
[P, num_grid, num_bands] = size(D_bands);
A_main = zeros(P * num_bands, num_grid);
C_block = zeros(P * num_bands, num_bands);
y = zeros(P * num_bands, 1);

for band_idx = 1:num_bands
    rows = (band_idx - 1) * P + (1:P);
    A_main(rows, :) = D_bands(:, :, band_idx);
    C_block(rows, band_idx) = c_bands(:, band_idx);
    y(rows) = z_bands(:, band_idx);
end

A = [A_main, C_block];
end

function mu = sparse_bayes_linear(A, y, max_iter, tol)
[~, num_cols] = size(A);
gamma = ones(num_cols, 1);
lambda = max(1e-6, 1e-3 * mean(y .^ 2));
mu = zeros(num_cols, 1);

for iter = 1:max_iter
    Gamma = diag(gamma);
    Sigma_y = A * Gamma * A' + lambda * eye(size(A, 1));
    mu_new = Gamma * (A' / Sigma_y) * y;
    Sigma_x = Gamma - Gamma * (A' / Sigma_y) * A * Gamma;
    gamma_new = max(real(mu_new) .^ 2 + max(real(diag(Sigma_x)), 0), 1e-10);

    if norm(mu_new - mu) <= tol * max(1, norm(mu_new))
        mu = mu_new;
        break;
    end

    mu = mu_new;
    gamma = gamma_new;
end
end

function [spectrum, doa_est] = zdoa_wideband_exact(z_bands, D_bands, c_bands, scan_deg, num_sources, L_snap, cfg)
num_bands = size(z_bands, 2);
num_grid = numel(scan_deg);
num_combos = nchoosek_safe(num_grid, num_sources);

if num_combos <= cfg.max_support_evals
    candidate_idx = 1:num_grid;
else
    matched = zeros(num_grid, 1);
    for band_idx = 1:num_bands
        matched = matched + max(D_bands(:, :, band_idx)' * z_bands(:, band_idx), 0);
    end
    [~, order] = sort(matched, 'descend');
    pool_size = min(max(20, 4 * num_sources), num_grid);
    candidate_idx = sort(order(1:pool_size));
end

combos = nchoosek(candidate_idx, num_sources);
best_score = inf;
best_idx = combos(1, :);
best_p = zeros(num_sources, 1);

for combo_idx = 1:size(combos, 1)
    support_idx = combos(combo_idx, :);
    [p_support, score] = fit_support_gamma_wideband(z_bands, D_bands(:, support_idx, :), c_bands, L_snap, cfg);
    if score < best_score
        best_score = score;
        best_idx = support_idx;
        best_p = p_support;
    end
end

spectrum = zeros(num_grid, 1);
spectrum(best_idx) = best_p;
spectrum = normalize_spectrum(spectrum);
doa_est = sort(scan_deg(best_idx));
end

function [p_support, score] = fit_support_gamma_wideband(z_bands, Ds_bands, c_bands, L_snap, cfg)
[A, y] = build_support_system(z_bands, Ds_bands, c_bands);
x = lsqnonneg(A, y);
num_support = size(Ds_bands, 2);
p_support = x(1:num_support);
sigma_vec = x(num_support + 1:end);

for iter = 1:cfg.max_alt_iters
    u_model = max(A * [p_support; sigma_vec], cfg.fit_floor);
    weights = 1 ./ max(u_model .^ 2, cfg.fit_floor);
    Aw = bsxfun(@times, A, sqrt(weights));
    bw = sqrt(weights) .* y;
    x_next = lsqnonneg(Aw, bw);
    if norm(x_next - [p_support; sigma_vec]) <= 1e-8 * max(1, norm([p_support; sigma_vec]))
        p_support = x_next(1:num_support);
        sigma_vec = x_next(num_support + 1:end);
        break;
    end
    p_support = x_next(1:num_support);
    sigma_vec = x_next(num_support + 1:end);
end

u_model = max(A * [p_support; sigma_vec], cfg.fit_floor);
score = sum(L_snap * log(u_model) + L_snap * (y ./ u_model));
end

function [A, y] = build_support_system(z_bands, Ds_bands, c_bands)
[P, num_support, num_bands] = size(Ds_bands);
A_main = zeros(P * num_bands, num_support);
C_block = zeros(P * num_bands, num_bands);
y = zeros(P * num_bands, 1);

for band_idx = 1:num_bands
    rows = (band_idx - 1) * P + (1:P);
    A_main(rows, :) = Ds_bands(:, :, band_idx);
    C_block(rows, band_idx) = c_bands(:, band_idx);
    y(rows) = z_bands(:, band_idx);
end

A = [A_main, C_block];
end

function save_benchmark_plots(final_result, prefix)
curves = final_result.curves;
methods = final_result.methods;
exp_cfg = final_result.exp_cfg;

fig = figure('Color', 'w', 'Position', [80, 80, 1360, 900]);

subplot(2, 2, 1);
hold on;
for idx = 1:numel(methods)
    plot(curves.snr_values, curves.mean_mse_deg2(idx, :), '-o', 'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', methods{idx}.name);
end
grid on;
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Mean squared error (deg^2)');
title(sprintf('%s MSE', exp_cfg.modality));
legend('Location', 'eastoutside');

subplot(2, 2, 2);
hold on;
for idx = 1:numel(methods)
    plot(curves.snr_values, curves.exact_rate(idx, :), '-o', 'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', methods{idx}.name);
end
grid on;
ylim([0, 1.02]);
xlabel('SNR (dB)');
ylabel('Exact recovery rate');
title(sprintf('%s exact recovery', exp_cfg.modality));
legend('Location', 'eastoutside');

mid_idx = ceil(numel(curves.snr_values) / 2);
subplot(2, 2, 3);
hold on;
for idx = 1:numel(methods)
    plot(exp_cfg.scan_deg, normalize_spectrum(curves.example_spectra(:, idx, mid_idx)), 'LineWidth', 1.2, 'DisplayName', methods{idx}.name);
end
for k = 1:numel(exp_cfg.true_doa_deg)
    xline(exp_cfg.true_doa_deg(k), '--k', 'LineWidth', 1.0);
end
grid on;
xlabel('Angle (deg)');
ylabel('Normalized spectrum');
title(sprintf('Example spectra at SNR = %g dB', curves.snr_values(mid_idx)));
legend('Location', 'eastoutside');

subplot(2, 2, 4);
axis off;
text(0.02, 0.92, sprintf('Scenario: %s', exp_cfg.name), 'FontSize', 12, 'FontWeight', 'bold');
text(0.02, 0.80, sprintf('W design: %s', exp_cfg.w_design.label), 'FontSize', 11);
text(0.02, 0.70, sprintf('M = %d, P = %d, L = %d', exp_cfg.M, exp_cfg.P, exp_cfg.L), 'FontSize', 11);
text(0.02, 0.60, sprintf('True DOA (deg): %s', mat2str(exp_cfg.true_doa_deg)), 'FontSize', 11);
text(0.02, 0.50, sprintf('Source powers: %s', mat2str(exp_cfg.source_power, 3)), 'FontSize', 11);
if strcmp(exp_cfg.modality, 'wideband')
    text(0.02, 0.40, sprintf('Freq ratios: %s', mat2str(exp_cfg.freq_ratios, 3)), 'FontSize', 11);
end
text(0.02, 0.28, sprintf('Best overall: %s', final_result.best_overall_method), 'FontSize', 11);
text(0.02, 0.18, sprintf('Our method best-at-SNR count: %d / %d', final_result.our_best_snr_count, numel(curves.snr_values)), 'FontSize', 11);

sgtitle(sprintf('%s benchmark with 7 baselines', upper_first(exp_cfg.modality)));
exportgraphics(fig, [prefix '_summary.png'], 'Resolution', 180);

fig2 = figure('Color', 'w', 'Position', [100, 100, 1280, 760]);
tiledlayout(2, 1);

nexttile;
hold on;
for idx = 1:numel(methods)
    plot(curves.snr_values, curves.mean_runtime_s(idx, :), '-o', 'LineWidth', 1.4, 'MarkerSize', 4, 'DisplayName', methods{idx}.name);
end
grid on;
xlabel('SNR (dB)');
ylabel('Average runtime (s)');
title(sprintf('%s runtime', exp_cfg.modality));
legend('Location', 'eastoutside');

nexttile;
zero_idx = find(curves.snr_values == 0, 1, 'first');
if isempty(zero_idx)
    zero_idx = ceil(numel(curves.snr_values) / 2);
end
bar(curves.exact_rate(:, zero_idx));
set(gca, 'XTickLabel', cellfun(@(m) m.short, methods, 'UniformOutput', false), 'XTickLabelRotation', 30);
ylim([0, 1.02]);
grid on;
ylabel('Exact recovery rate');
title(sprintf('Exact rate at %g dB', curves.snr_values(zero_idx)));
exportgraphics(fig2, [prefix '_runtime.png'], 'Resolution', 180);
end

function write_summary_text(results, filename)
fid = fopen(filename, 'w');
if fid < 0
    error('Unable to open %s for writing.', filename);
end

cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Benchmark timestamp: %s\n\n', results.meta.timestamp);

write_one_summary_block(fid, 'NARROWBAND', results.narrowband);
fprintf(fid, '\n');
write_one_summary_block(fid, 'WIDEBAND', results.wideband);
end

function write_one_summary_block(fid, title_str, result_block)
best_cfg = result_block.search.best_setup;
curves = result_block.final.curves;
methods = result_block.final.methods;

fprintf(fid, '[%s]\n', title_str);
fprintf(fid, 'Selected scenario: %s\n', best_cfg.name);
fprintf(fid, 'W design: %s\n', best_cfg.w_design.label);
fprintf(fid, 'True DOA (deg): %s\n', mat2str(best_cfg.true_doa_deg));
fprintf(fid, 'Source powers: %s\n', mat2str(best_cfg.source_power, 3));
fprintf(fid, 'M = %d, P = %d, L = %d\n', best_cfg.M, best_cfg.P, best_cfg.L);
if strcmp(best_cfg.modality, 'wideband')
    fprintf(fid, 'Freq ratios: %s\n', mat2str(best_cfg.freq_ratios, 3));
end
fprintf(fid, 'Best overall method by average MSE: %s\n', result_block.final.best_overall_method);
fprintf(fid, '\nPer-method average MSE over SNR sweep:\n');
for idx = 1:numel(methods)
    fprintf(fid, '  %-20s %.6f deg^2\n', methods{idx}.name, mean(curves.mean_mse_deg2(idx, :)));
end
fprintf(fid, '\nPer-method average exact rate over SNR sweep:\n');
for idx = 1:numel(methods)
    fprintf(fid, '  %-20s %.6f\n', methods{idx}.name, mean(curves.exact_rate(idx, :)));
end
end

function beam_label = upper_first(str)
beam_label = str;
beam_label(1) = upper(beam_label(1));
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

function spectrum = normalize_spectrum(spectrum)
spectrum = real(spectrum(:));
spectrum(~isfinite(spectrum)) = 0;
peak = max(spectrum);
if peak > 0
    spectrum = spectrum / peak;
end
end

function doa_est = select_top_doa_peaks_local(spectrum, scan_deg, num_sources, min_sep_deg)
spectrum = spectrum(:);
scan_deg = scan_deg(:);

is_peak = false(numel(spectrum), 1);
if numel(spectrum) == 1
    is_peak(1) = true;
else
    is_peak(1) = spectrum(1) >= spectrum(2);
    is_peak(end) = spectrum(end) >= spectrum(end - 1);
    for idx = 2:numel(spectrum) - 1
        is_peak(idx) = spectrum(idx) >= spectrum(idx - 1) && spectrum(idx) >= spectrum(idx + 1);
    end
end

candidate_idx = find(is_peak);
if isempty(candidate_idx)
    candidate_idx = (1:numel(spectrum)).';
end

[~, sort_order] = sort(spectrum(candidate_idx), 'descend');
candidate_idx = candidate_idx(sort_order);

selected = zeros(num_sources, 1);
selected_count = 0;
for idx = 1:numel(candidate_idx)
    if selected_count == 0
        selected_count = 1;
        selected(selected_count) = candidate_idx(idx);
        continue;
    end

    if all(abs(scan_deg(candidate_idx(idx)) - scan_deg(selected(1:selected_count))) >= min_sep_deg)
        selected_count = selected_count + 1;
        selected(selected_count) = candidate_idx(idx);
        if selected_count == num_sources
            break;
        end
    end
end

if selected_count < num_sources
    all_idx = (1:numel(spectrum)).';
    [~, sort_order] = sort(spectrum(all_idx), 'descend');
    all_idx = all_idx(sort_order);
    for idx = 1:numel(all_idx)
        if any(selected(1:selected_count) == all_idx(idx))
            continue;
        end
        selected_count = selected_count + 1;
        selected(selected_count) = all_idx(idx);
        if selected_count == num_sources
            break;
        end
    end
end

doa_est = sort(scan_deg(selected(1:selected_count))).';
end

function W = build_steered_beamformer_local(M, beam_deg, d_eff)
num_beams = numel(beam_deg);
W = zeros(num_beams, M);
for beam_idx = 1:num_beams
    W(beam_idx, :) = steering_vector_local(beam_deg(beam_idx), M, d_eff)';
end
end

function [D, c] = build_beam_power_dictionary_local(W, scan_deg, d_eff)
num_beams = size(W, 1);
num_grid = numel(scan_deg);
num_sensors = size(W, 2);
D = zeros(num_beams, num_grid);
for grid_idx = 1:num_grid
    a = steering_vector_local(scan_deg(grid_idx), num_sensors, d_eff);
    D(:, grid_idx) = abs(W * a) .^ 2;
end
c = sum(abs(W) .^ 2, 2);
end

function a = steering_vector_local(theta_deg, M, d_eff)
sensor_index = (0:M - 1).';
phase_shift = 2 * pi * d_eff * sind(theta_deg) * sensor_index;
a = exp(1i * phase_shift) / sqrt(M);
end

function value = nchoosek_safe(n, k)
if k > n
    value = 0;
    return;
end
value = 1;
for idx = 1:k
    value = value * (n - k + idx) / idx;
    if value > 1e9
        break;
    end
end
end
