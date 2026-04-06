%% ============ LNCA算法测试脚本（修复版）=============
%% 基于论文黄-2022 P84-91的l-范数循环算法

clear; clc; close all;

%% =============参数设置=============
N = 256;               % 总子载波数
Nc = 64;               % 通信子载波数
Nr = N - Nc;           % 雷达子载波数
modSC = 'QPSK';        % 调制方式
e = 0.1;               % 通信能量占总能量比例
monte = 10;           % Monte Carlo次数
l_norm_order = 100;     % l-范数阶数（推荐10~100）
iteration = 100;        % LNCA迭代次数

rng(1)  % 固定随机种子

% ===== BER Eb/N0 配置 =====
ebno_list = 0:1:18;       % Eb/N0 范围（论文标准）
nEbNo = length(ebno_list);
bits_per_symbol = 2;      % QPSK = 2 bit


%% =============波形配置=============
waveforms = {
    'OFDM原始', '原始OFDM';
    'LNCA优化', 'l-范数循环算法';
    'OCDM','OCDM';
};

nWave = length(waveforms);
Progress = waitbar(0, 'Progress...');

%% =============预分配结果存储============
papr_results = cell(1, nWave);      % PAPR结果
aisl_results = cell(1, nWave);      % AISL结果（可选）
ber_all = cell(nWave, nEbNo);  % ber_all{波形, Eb/N0点}

for w = 1:nWave
    papr_results{w} = zeros(1, monte);
    aisl_results{w} = zeros(1, monte);
end

%% =============辅助函数=============
PAPR_fun = @(x) max(abs(x(:)).^2) / mean(abs(x(:)).^2);
FH = @(x) ifft(x)*sqrt(N);
F = @(x) fft(x)*sqrt(N);
function DFnT_matrix=DFnT(P)
    m=(0:P-1)'*ones(1,P);
    n=ones(P,1)*(0:P-1);

    fresnel_kernel=exp(1j*pi/P*(m-n).^2);
    normalization=(1/sqrt(P))*exp((pi/4)*1j);
    DFnT_matrix=normalization*fresnel_kernel;
end

%% =============Monte Carlo循环============
for m = 1:monte

    % 进度条
    waitbar(m/monte, Progress, sprintf('Progress: %d/%d (%.1f%%)', m, monte, m/monte*100));

    % ===== 1. 随机子载波划分 =====
    rand_idx = randperm(N);
    idx_comm = sort(rand_idx(1:Nc));
    mask_comm = false(N, 1);
    mask_comm(idx_comm) = true;
    idx_radar = find(~mask_comm);

    % ===== 2. 生成通信码（恒模）=====
    switch upper(modSC)
        case 'BPSK'
            bits = randi([0 1], Nc, 1);
            zc_unit = exp(1j * pi * bits);
        case 'QPSK'
            bits = randi([0 1], 2*Nc, 1);
            zc_unit = (1/sqrt(2)) * ((2*(bits(1:2:end)==1)-1) + 1j*(2*(bits(2:2:end)==1)-1));
        otherwise
            error('仅支持 BPSK/QPSK');
    end

    % 通信码能量归一化
    zc = sqrt(e) * zc_unit / sqrt(Nc);

    % ===== 3. 雷达码（自由变量 + 归一化）=====
    Er = 1 - e;
    if Er <= 0
        error('e 设置过大，导致无可用雷达能量');
    end

    zr = randn(Nr, 1) + 1j * randn(Nr, 1);   % 初始随机
    zr = zr / norm(zr) * sqrt(Er);           % 能量归一化

    % ===== 4. 合成频域信号 =====
    z = zeros(N, 1);
    z(idx_comm) = zc;
    z(idx_radar) = zr;

    % ===== 5. 波形生成与优化 =====
    for w = 1:nWave
        switch w
            case 1  % 原始OFDM
                % 4倍过采样IFFT
                z_padded = [z; zeros(3*N, 1)];
                x = ifft(z_padded) * sqrt(4*N);
                time_signal = x(1:4*N);  % 取4N点

            case 2  % LNCA优化
                % 调用修复后的LNCA算法
                z_opt = l_norm_fixed(z, idx_comm, idx_radar, l_norm_order, iteration);

                % 4倍过采样IFFT生成时域信号
                z_padded = [z_opt; zeros(3*N, 1)];
                x = ifft(z_padded) * sqrt(4*N);
                time_signal = x(1:4*N);

            case 3 % OCDM
                Psi4 = DFnT(4*N);  % 4倍过采样 DFnT
                z_padded = zeros(4*N, 1);
                z_padded(1:N) = z;
                time_signal = Psi4' * z_padded * sqrt(4*N);
            case 4
                time_signal = optimize_waveform_cvx(z, idx_comm, idx_radar);
            case 5
                time_signal = FH(optimize_waveform_freq(z, idx_comm, idx_radar));

        end

        % 计算PAPR
        papr_results{w}(m) = PAPR_fun(time_signal);
    end

end

%% =============CCDF绘制=============

thresholds = 0:0.1:12;

figure;
colors = {'b-', 'r-', 'g-', 'm-', 'c-'};
line_widths = [1.5, 2.0];

for w = 1:nWave
    papr_dB = 10 * log10(papr_results{w});
    ccdf = arrayfun(@(t) mean(papr_dB > t), thresholds);
    semilogy(thresholds, ccdf, colors{w}, 'LineWidth', line_widths(min(w, length(line_widths))));
    hold on;
end


grid on;
xlabel('PAPR (dB)', 'FontSize', 12);
ylabel('CCDF', 'FontSize', 12);
title(sprintf('CCDF of PAPR (N=%d, Nc=%d, l=%d, Iter=%d, MC=%d)', N, Nc, l_norm_order, iteration, monte), 'FontSize', 12);
legend(waveforms(:, 1), 'Location', 'northeast', 'FontSize', 10);

% 标记关键点（CCDF=1e-4处的PAPR）
for w = 1:nWave
    papr_dB = 10 * log10(papr_results{w});
    target_ccdf = 1e-4;
    [~, idx] = min(abs(arrayfun(@(t) mean(papr_dB > t), thresholds) - target_ccdf));
    if idx > 0 && idx <= length(thresholds)
        fprintf('%s @ CCDF=1e-4: PAPR ≈ %.2f dB\n', waveforms{w,1}, thresholds(idx));
    end
end

saveas(gcf, sprintf('CCDF_N%d_Nc%d_e%.1f.png', N, Nc, e));
savefig(gcf, sprintf('CCDF_N%d_Nc%d_e%.1f.fig', N, Nc, e));
hold off;

fprintf('Done!\n');