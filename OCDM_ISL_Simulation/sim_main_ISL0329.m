%% ============参数=============
N = 256;
Nc = 64;
Nr = N - Nc;
modSC = 'QPSK';
monte = 1;        % Monte Carlo次数
e = 0.3;              % 通信总能量
thresholds = 0:0.2:8;
rng(32);

% ===== 波形配置 =====
% 格式: {名称, 变换函数句柄}
% 变换函数输入频域信号z，输出时域信号x
waveforms = {
    'OFDM',         @(z) ifft(z) * sqrt(N);
    'OCDM',         @(z) DFnT(N)' * z;
    'optimize_waveform_cvx',@(z) ifft(optimize_waveform_cvx(z, idx_comm, idx_radar)) * sqrt(N);
    'optimize_waveform_freq', @(z) ifft(optimize_waveform_freq(z, idx_comm, idx_radar)) * sqrt(N);
    'l_norm_fixed' @(z) ifft(l_norm_fixed(z, idx_comm, idx_radar, 20 ,50)) * sqrt(N);
};

nWave = length(waveforms);
Progress = waitbar(0, 'Progress...');

% ===== 预分配结果存储 =====
papr_all = cell(1, nWave);      % PAPR结果

for w = 1:nWave
    papr_all{w} = zeros(1, monte);
end

% 保存最后一次的频域信号z（用于画自相关）
z_final = zeros(N, 1);

%% =============函数============
PAPR_fun = @(x) max(abs(x(:)).^2) / mean(abs(x(:)).^2);

%DFnT matrix
function DFnT_matrix=DFnT(P)
    m=(0:P-1)'*ones(1,P);
    n=ones(P,1)*(0:P-1);

    fresnel_kernel=exp(1j*pi/P*(m-n).^2);
    normalization=(1/sqrt(P))*exp((pi/4)*1j);
    DFnT_matrix=normalization*fresnel_kernel;
end


%% ============monte carlo循环===========
for m = 1:monte

    %进度条
    waitbar(m/monte, Progress, sprintf('Progress: %d/%d (%.1f%%)', m, monte, m/monte*100));

    % ===== 1. 随机子载波划分 =====
    rand_idx = randperm(N);
    idx_comm = sort(rand_idx(1:Nc));
    mask_comm = false(N,1);
    mask_comm(idx_comm) = true;
    idx_radar = find(~mask_comm);

    % ===== 2. 通信码（恒模 = e）=====
    switch upper(modSC)
        case 'BPSK'
            zc_unit = exp(1j*pi*randi([0 1], Nc, 1));
        case 'QPSK'
            bits = randi([0 1], 2*Nc, 1);
            zc_unit = (1/sqrt(2))*((2*(bits(1:2:end)==1)-1) + 1j*(2*(bits(2:2:end)==1)-1));
        otherwise
            error('仅支持 BPSK/QPSK');
    end

    zc = sqrt(e) * zc_unit / sqrt(Nc);

    Er = 1 - e;

    if Er <= 0
        error('e 设置过大，导致无可用雷达能量');
    end

    % ===== 4. 雷达码（自由变量 + 归一化）=====
    zr = randn(Nr,1) + 1j*randn(Nr,1);   % 初始随机
    zr = zr / norm(zr) * sqrt(Er);       % 能量归一化

    % ===== 5. 合成频域信号 =====
    z = zeros(N,1);
    z(idx_comm)  = zc;
    z(idx_radar) = zr;

    % ===== 6. 遍历各波形 =====
    for w = 1:nWave
        transform_func = waveforms{w, 2};
        switch w
            case 1
                time_signal = transform_func(z);  % 统一调用变换函数
            case 2
                time_signal = transform_func(z);
            case 3
                time_signal = transform_func(z);
            case 4
                time_signal = transform_func(z);
            case 5
                time_signal = transform_func(z);
        end
        papr_all{w}(m) = PAPR_fun(time_signal);
    end

    % 保存最后一次的频域信号z（用于画单次自相关）
    if m == monte
        z_final = z;
    end

end


%% =============CCDF=============
figure;
colors = {'b-', 'r-', 'g-', 'm-', 'c-'};  % 颜色配置

for w = 1:nWave
    papr_dB = 10*log10(papr_all{w});
    ccdf = arrayfun(@(t) mean(papr_dB > t), thresholds);
    semilogy(thresholds, ccdf, colors{w}, 'LineWidth', 1.5);
    hold on;
end
hold off;
grid on;
xlabel('PAPR (dB)');
ylabel('CCDF');
title('CCDF of PAPR');
legend(waveforms(:,1), 'Location', 'northeast');


%% ==============自相关==============
% 使用最后一次Monte Carlo的波形来计算自相关
figure;
lags = -(N-1):(N-1);

for w = 1:nWave
    transform_func = waveforms{w, 2};
    x = transform_func(z_final);    % 用最后一次的z变换
    
    r = xcorr(x);
    r0 = r(N);
    r_dB = 20*log10(max(abs(r)/max(abs(r0),1e-12),1e-12));
    plot(lags, r_dB, colors{w}, 'LineWidth', 2);
    hold on;
end
hold off;
grid on;
xlabel('Time Lag');
ylabel('Autocorrelation Level (dB)');
title('自相关特性（单次波形）');
ylim([-80 0]);
xlim([-(N-1) N-1]);
legend(waveforms(:,1), 'Location', 'northeast');

close(Progress);
disp('Done!');
