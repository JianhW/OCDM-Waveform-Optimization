%% OCDM / OFDM CCDF simulation
% Compare baseline and optimization methods under a unified CCDF setup.

clear; clc; close all;

%% Parameters
N = 256;                   % Total number of subcarriers
Nc = 32;                   % Number of communication subcarriers
Nr = N - Nc;               % Number of radar subcarriers
modSC = 'QPSK';            % Communication modulation
e = 0.1;                   % Communication energy ratio
monte = 10;                % Monte Carlo trials
l_norm_order = 100;        % l-norm order
iteration = 1000;          % Main iteration number
os_factor = 4;             % Oversampling factor for CCDF evaluation

% Parameters for the proposed OCDM algorithm in Mine.m
mine_lambda = 0.0016;
mine_warm_start_iteration = 500;
mine_pisl_guard = 1e-7;
mine_lnorm_bound = 'current';
mine_eig_mode = 'upper';
mine_max_backtracking = 14;
mine_lambda_shrink = 0.5;

rng(1);


%% Waveform list
waveforms = {
    'OFDM',      'Original OFDM';
    'OCDM',      'Original OCDM';
    'Mine',      'Joint PISL-lNorm OCDM';
    'Varshney',  'Varshney';
    'Wang',      'Wang';
};

nWave = size(waveforms, 1);
Progress = waitbar(0, 'Progress...');

%% Storage
papr_results = cell(1, nWave);
aisl_results = cell(1, nWave); 

for w = 1:nWave
    papr_results{w} = zeros(1, monte);
    aisl_results{w} = zeros(1, monte);
end

%% =============辅助函数=============
PAPR_fun = @(x) max(abs(x(:)).^2) / mean(abs(x(:)).^2);
FH = @(x) ifft(x)*sqrt(N);
F = @(x) fft(x)*sqrt(N);

function DFnT_matrix = DFnT(P)
    m = (0:P-1)' * ones(1, P);
    n = ones(P, 1) * (0:P-1);
    fresnel_kernel = exp(1j * pi / P * (m - n).^2);
    normalization = (1 / sqrt(P)) * exp((pi / 4) * 1j);
    DFnT_matrix = normalization * fresnel_kernel;
end

Psi_os = DFnT(os_factor * N);

%% Monte Carlo loop
for m = 1:monte

    waitbar(m / monte, Progress, sprintf('Progress: %d/%d (%.1f%%)', m, monte, 100 * m / monte));

    % ===== 1. Random subcarrier allocation =====
    rand_idx = randperm(N);
    idx_comm = sort(rand_idx(1:Nc)).';
    mask_comm = false(N, 1);
    mask_comm(idx_comm) = true;
    idx_radar = find(~mask_comm);

    % ===== 2. Communication symbols =====
    switch upper(modSC)
        case 'BPSK'
            bits = randi([0 1], Nc, 1);
            zc_unit = exp(1j * pi * bits);
        case 'QPSK'
            bits = randi([0 1], 2 * Nc, 1);
            zc_unit = (1 / sqrt(2)) * ...
                ((2 * bits(1:2:end) - 1) + 1j * (2 * bits(2:2:end) - 1));
        otherwise
            error('Only BPSK/QPSK are supported.');
    end

    zc = sqrt(e) * zc_unit / sqrt(Nc);

    % ===== 3. Radar symbols =====
    Er = 1 - e;
    if Er <= 0
        error('Parameter e is too large and leaves no radar energy.');
    end

    zr = randn(Nr, 1) + 1j * randn(Nr, 1);
    zr = zr / norm(zr) * sqrt(Er);

    % ===== 4. Composite frequency-domain vector =====
    z = zeros(N, 1);
    z(idx_comm) = zc;
    z(idx_radar) = zr;

    % ===== 5. Waveform generation / optimization =====
    for w = 1:nWave
        switch w
            case 1  % Original OFDM
                z_padded = [z; zeros((os_factor - 1) * N, 1)];
                time_signal = ifft(z_padded) * sqrt(os_factor * N);

            case 2  % Original OCDM
                z_padded = zeros(os_factor * N, 1);
                z_padded(1:N) = z;
                time_signal = Psi_os' * z_padded * sqrt(os_factor * N);

            case 3  % Proposed algorithm
                z_mine = Mine(z, idx_comm, idx_radar, ...
                    'maxIt', iteration, ...
                    'minIt', 200, ...
                    'lNormOrder', l_norm_order, ...
                    'lambda', mine_lambda, ...
                    'osFactor', os_factor, ...
                    'epsStop', mine_pisl_guard, ...
                    'lnormBound', mine_lnorm_bound, ...
                    'etaMode', mine_eig_mode, ...
                    'maxBacktracking', mine_max_backtracking, ...
                    'backtrackingShrink', mine_lambda_shrink, ...
                    'UpdateMode', 'hybrid', ...
                    'PaprWeight', 0.85);
                z_padded = zeros(os_factor * N, 1);
                z_padded(1:N) = z_mine;
                time_signal = Psi_os' * z_padded * sqrt(os_factor * N);

            case 4  % Varshney
                time_signal = Varshney(z, idx_comm, idx_radar);

            case 5  % Wang
                time_signal = FH(Wang(z, idx_comm, idx_radar));
        end

        papr_results{w}(m) = PAPR_fun(time_signal);
    end
end

%% CCDF plotting
if isgraphics(Progress)
    close(Progress);
end

thresholds = 0:0.1:12;
ccdf_curves = cell(1, nWave);

for w = 1:nWave
    papr_dB = 10 * log10(papr_results{w});
    ccdf_curves{w} = arrayfun(@(t) mean(papr_dB > t), thresholds);
end

plot_colors = [
    0.0000, 0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.4660, 0.6740, 0.1880;
    0.4940, 0.1840, 0.5560;
    0.3010, 0.7450, 0.9330
];
line_styles = {'-', '--', '-.', ':', '-'};
marker_styles = {'o', 's', 'd', '^', 'v'};
line_width = 3.0;
marker_size = 7.5;
n_markers = 9;

legend_handles = gobjects(1, nWave);

fig = figure('Color', 'w', 'Position', [120, 120, 760, 560]);
ax = axes(fig);
hold(ax, 'on');
box(ax, 'on');
set(ax, 'FontName', 'Times New Roman', ...
        'FontSize', 13, ...
        'LineWidth', 1.2, ...
        'TickDir', 'in', ...
        'TickLength', [0.018, 0.018], ...
        'YScale', 'log');

for w = 1:nWave
    ccdf_plot = ccdf_curves{w};
    ccdf_plot(ccdf_plot == 0) = NaN;

    valid_idx = find(~isnan(ccdf_plot));
    
    if ~isempty(valid_idx)

        x_valid = thresholds(valid_idx);
    
        y_valid = log10(ccdf_plot(valid_idx));
   
        dx = diff(x_valid);
        dy = diff(y_valid);
    
        seg_len = sqrt(dx.^2 + dy.^2);
    
        cum_len = [0 cumsum(seg_len)];

        target_len = linspace(0, cum_len(end), ...
            min(n_markers, length(valid_idx)));

        marker_local_idx = arrayfun(@(s) ...
            find(abs(cum_len - s) == min(abs(cum_len - s)), 1), ...
            target_len);

        marker_idx = valid_idx(unique(marker_local_idx));
    
    else
        marker_idx = [];
    
    end

    legend_handles(w) = semilogy(ax, thresholds, ccdf_plot, ...
        'LineStyle', line_styles{mod(w - 1, length(line_styles)) + 1}, ...
        'Color', plot_colors(mod(w - 1, size(plot_colors, 1)) + 1, :), ...
        'LineWidth', line_width, ...
        'Marker', marker_styles{mod(w - 1, length(marker_styles)) + 1}, ...
        'MarkerIndices', marker_idx, ...
        'MarkerSize', marker_size, ...
        'MarkerFaceColor', 'w');
end

grid(ax, 'on');
ax.GridAlpha = 0.18;
ax.MinorGridAlpha = 0.10;
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.Layer = 'top';
xlim(ax, [thresholds(1), thresholds(end)]);

valid_curves = ccdf_curves(cellfun(@(c) any(c > 0), ccdf_curves));
if ~isempty(valid_curves)
    min_positive_ccdf = min(cellfun(@(c) min(c(c > 0)), valid_curves));
    ylim(ax, [10^(floor(log10(min_positive_ccdf))), 1]);
end

xlabel(ax, 'PAPR (dB)', 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax, 'CCDF', 'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold');
title(ax, sprintf('CCDF of PAPR (N = %d, N_c = %d, l = %d, Iter = %d, MC = %d)', ...
    N, Nc, l_norm_order, iteration, monte), ...
    'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');

lgd = legend(ax, legend_handles, waveforms(:, 1), ...
    'Location', 'southwest', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 11);
lgd.Box = 'off';

target_ccdf = 1e-4;
for w = 1:nWave
    ccdf_current = ccdf_curves{w};
    [~, idx] = min(abs(ccdf_current - target_ccdf));
    if ~isempty(idx) && idx >= 1 && idx <= length(thresholds)
        fprintf('%s @ nearest CCDF to 1e-4: PAPR = %.2f dB, CCDF = %.4g\n', ...
            waveforms{w, 1}, thresholds(idx), ccdf_current(idx));
    end
end

print(fig, sprintf('CCDF_N%d_Nc%d_e%.1f_lambda%.3g.png', N, Nc, e, mine_lambda), '-dpng', '-r600');
savefig(fig, sprintf('CCDF_main_first1.fig'));
hold(ax, 'off');

fprintf('Done!\n');
