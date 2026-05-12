%% ============参数=============
N = 256;
Nc = 64;
Nr = N - Nc;
modSC = 'QPSK';
monte = 1;        % Monte Carlo次数
e = 0.2;              % 通信总能量
rng(32);
pisl_iteration = 500;

% ===== 波形配置 =====
% 格式: {名称, 变换函数句柄}
% 变换函数输入频域信号z，输出时域信号x
waveforms = {
    'OFDM',         @(z) ifft(z) * sqrt(N);
    'OCDM',         @(z) DFnT(N)' * z;
    'optimize_waveform_cvx',@(z) ifft(optimize_waveform_cvx(z, idx_comm, idx_radar)) * sqrt(N);
    'optimize_waveform_freq', @(z) ifft(optimize_waveform_freq(z, idx_comm, idx_radar)) * sqrt(N);
    'OCDM-PISL-MM' @(z) DFnT(N)' * OCDM_PISL_MM(z, idx_comm, idx_radar, pisl_iteration, ...
        'Transform', 'inverse', 'TotalNorm', 1);
};

nWave = length(waveforms);
Progress = waitbar(0, 'Progress...');

% ===== 预分配结果存储 =====
pisl_all = cell(1, nWave);      % PISL结果

for w = 1:nWave
    pisl_all{w} = zeros(1, monte);
end

% 保存最后一次的频域信号z（用于画自相关）
z_final = zeros(N, 1);

%% =============函数============
PISL_fun = @(x) calc_PISL(x);


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
        pisl_all{w}(m) = PISL_fun(time_signal);
    end

    % 保存最后一次的频域信号z（用于画单次自相关）
    if m == monte
        z_final = z;
    end

end


%% =============PISL对比=============
figure;
colors = {'b-', 'r-', 'g-', 'm-', 'c-', 'k-'};  % 颜色配置

pisl_mean = zeros(1, nWave);
for w = 1:nWave
    pisl_mean(w) = mean(pisl_all{w});
end
bar(pisl_mean);
grid on;
set(gca, 'XTickLabel', waveforms(:,1));
ylabel('Mean PISL');
title(sprintf('PISL comparison (N=%d, Nc=%d, MC=%d)', N, Nc, monte));
legend(waveforms(:,1), 'Location', 'northeast');

for w = 1:nWave
    fprintf('%s: Mean PISL = %.6e\n', waveforms{w,1}, pisl_mean(w));
end


%% ==============自相关==============
% 使用最后一次Monte Carlo的波形来计算自相关
figure;
lags = -(N-1):(N-1);

for w = 1:nWave
    transform_func = waveforms{w, 2};
    x = transform_func(z_final);    % 用最后一次的z变换
    
    r = periodic_autocorr(x);
    r_bilateral = [conj(flipud(r(2:end))); r];
    r0 = r_bilateral(N);
    r_dB = 20*log10(max(abs(r_bilateral)/max(abs(r0),1e-12),1e-12));
    plot(lags, r_dB, colors{w}, 'LineWidth', 2);
    hold on;
end
hold off;
grid on;
xlabel('Time Lag');
ylabel('Periodic Autocorrelation Level (dB)');
title('周期自相关特性（单次波形）');
ylim([-80 0]);
xlim([-(N-1) N-1]);
legend(waveforms(:,1), 'Location', 'northeast');

close(Progress);
disp('Done!');

function pisl = calc_PISL(x)
    x = x(:);
    r = periodic_autocorr(x);
    r(1) = 0;
    pisl = 2 * sum(abs(r).^2);
    pisl = real(pisl);
end

function r = periodic_autocorr(x)
    x = x(:);
    r = ifft(abs(fft(x)).^2);
end

function DFnT_matrix=DFnT(P)
    m=(0:P-1)'*ones(1,P);
    n=ones(P,1)*(0:P-1);

    fresnel_kernel=exp(1j*pi/P*(m-n).^2);
    normalization=(1/sqrt(P))*exp(1j*pi/4);
    DFnT_matrix=normalization*fresnel_kernel;
end
