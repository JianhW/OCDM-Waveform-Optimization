%% ============参数=============
N = 256;
Nc = 16;
Nr = N - Nc;
modSC = 'QPSK';
monte = 10000;        % Monte Carlo次数
e = 0.5;              % 通信总能量
thresholds = 0:0.2:8;
alogrithms = {'OFDM', 'OCDM', 'OFDM New-ICF', 'OFDM Varshney'};

papr_OFDM_all = zeros(1, monte);  %记录papr
r_OFDM_acc = zeros(2*N-1,1);      %记录ISL

papr_OCDM_all = zeros(1, monte);  %记录papr
r_OCDM_acc = zeros(2*N-1,1);      %记录ISL

%% =============函数============
PAPR_fun = @(x) max(abs(x(:)).^2) / mean(abs(x(:)).^2);
FH = @(x) ifft(x)*sqrt(N);

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

    zc = sqrt(e) * zc_unit / sqrt(Nc);   % 恒模 = e

    % ===== 3. 剩余能量计算 =====
    Ec = sum(abs(zc).^2);     % 通信能量
    Er = 1 - Ec;              % 雷达可用能量

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

    x_OFDM = FH(z);
    x_OCDM = DFnT(N)'*z;

    

    % ===== 8. PAPR =====
    papr_OFDM_all(m) = PAPR_fun(x_OFDM);
    papr_OCDM_all(m) = PAPR_fun(x_OCDM);

    r_OFDM_tmp=xcorr(x_OFDM);
    r_OCDM_tmp=xcorr(x_OCDM);

    r_OFDM_acc=r_OFDM_acc+r_OFDM_tmp;
    r_OCDM_acc=r_OCDM_acc+r_OCDM_tmp;

end


%% =============CCDF=============
papr_OFDM_dB = 10*log10(papr_OFDM_all);
papr_OCDM_dB = 10*log10(papr_OCDM_all);


ccdf_OFDM = arrayfun(@(t) mean(papr_OFDM_dB > t), thresholds);
ccdf_OCDM = arrayfun(@(t) mean(papr_OCDM_dB > t), thresholds);

figure;
semilogy(thresholds, ccdf_OFDM,'LineWidth',1.5);
hold on;
semilogy(thresholds, ccdf_OCDM,'LineWidth',1.5);
hold off;
grid on;
xlabel('PAPR (dB)');
ylabel('CCDF');
title('CCDF of PAPR');
legend('OFDM','OCDM');

%% ==============xcorr==============
lags=-(N-1):(N-1);
r_avg = r_OFDM_acc / monte;
r0 = r_avg(N);
r_dB = 20*log10(max(abs(r_avg)/max(abs(r0),1e-12),1e-12));
%绘图
figure;
plot(lags, r_dB, 'LineWidth',2);
grid on;
xlabel('Time Lag');
ylabel('Autocorrelation Level (dB)');
title('自相关特性');
ylim([-80 0]);   % 可根据旁瓣大小调整
xlim([-(N-1) N-1]);
legend('OFDM');