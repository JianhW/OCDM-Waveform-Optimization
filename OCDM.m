clear; clc; close all;
N=128; %Number of Subcarrier
L=8; %Channel Length
Block_Num=8; %Block Number
C=16; %Len Cyclic Prefix
P=N+C;
monte=2;
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
modulation = [4];  %调制方式
numMod = length(modulation);

methods = {'OCDM_raw','OFDM_raw'};  % 一次可运行多个算法
numMethod = numel(methods);
methodCfgs = cell(1, numMethod);
for methodIdx = 1:numMethod
    methodCfgs{methodIdx} = get_method_cfg(methods{methodIdx});
end

equalizerList = 2;             % 1 = ZF, 2 = MMSE；这里只保留 MMSE
eqNames = {'MMSE'};
numEq = numel(equalizerList);

EbN0_dB = 0:4:40;
numSNR = numel(EbN0_dB);

total = zeros(numMethod, numSNR, numEq, numMod);  % 方法 × SNR × 均衡 × 调制
ratio = zeros(numMethod, numSNR, numEq, numMod);

for s = 1:numSNR
    dB = EbN0_dB(s);
    disp(dB);
    SNR = 10^(dB/10);

    for methodIdx = 1:numMethod
        methodCfg = methodCfgs{methodIdx};

        for eqIdx = 1:numEq
            Equal = equalizerList(eqIdx);

            for m = 1:numMod
                U = modulation(m);

                for loop = 1:monte

                    ber = process(methodCfg, U, Block_Num, N, C, L, Equal, SNR);
                    total(methodIdx, s, eqIdx, m) = total(methodIdx, s, eqIdx, m) + ber;
                    
                end
            end
        end
    end
end

total = total / monte;

figure();
ax = axes('FontName','Times New Roman', 'FontSize', 13, 'LineWidth', 1.1);
box(ax,'on'); hold(ax,'on'); grid(ax,'on');
set(ax, 'YScale', 'log', ...
    'GridColor', [0.82 0.82 0.82], ...
    'GridAlpha', 0.35, ...
    'MinorGridAlpha', 0.18, ...
    'XMinorGrid', 'off', ...
    'YMinorGrid', 'on', ...
    'Layer', 'top');

styleMap = get_plot_style_map();
legendNames = {};
curveIdx = 1;

for methodIdx = 1:numMethod
    methodCfg = methodCfgs{methodIdx};
    for eqIdx = 1:numEq
        for m = 1:numMod
            U = modulation(m);
            y = squeeze(total(methodIdx,:,eqIdx,m));
            styleCfg = styleMap{mod(curveIdx-1, numel(styleMap))+1};

            semilogy(ax, EbN0_dB, y, ...
                'Color', styleCfg.color, ...
                'LineStyle', styleCfg.lineStyle, ...
                'Marker', styleCfg.marker, ...
                'LineWidth', 2.3, ...
                'MarkerSize', 8.5, ...
                'MarkerIndices', 1:2:numSNR, ...
                'DisplayName', sprintf('%s-%s, %d-QAM', methodCfg.displayName, eqNames{eqIdx}, U));

            legendNames{end+1} = sprintf('%s-%s, %d-QAM', methodCfg.displayName, eqNames{eqIdx}, U);
            curveIdx = curveIdx + 1;
        end
    end
end

xlim(ax, [EbN0_dB(1), EbN0_dB(end)]);
xticks(ax, EbN0_dB);
ylim(ax, [1e-5 1]);
xlabel(ax, 'E_b/N_0 (dB)', 'FontName', 'Times New Roman', 'FontSize', 14);
ylabel(ax, 'Bit Error Rate', 'FontName', 'Times New Roman', 'FontSize', 14);
legend(ax, legendNames, 'Location', 'southwest', 'Interpreter', 'none', ...
    'FontName', 'Times New Roman', 'FontSize', 11, 'Box', 'off');

function styleMap = get_plot_style_map()
    styleMap = {
        struct('color', [0.00, 0.45, 0.74], 'lineStyle', '-',  'marker', 'o'), ...
        struct('color', [0.85, 0.33, 0.10], 'lineStyle', '-', 'marker', 's'), ...
        struct('color', [0.93, 0.69, 0.13], 'lineStyle', '-', 'marker', '^'), ...
        struct('color', [0.49, 0.18, 0.56], 'lineStyle', '-',  'marker', 'd'), ...
        struct('color', [0.47, 0.67, 0.19], 'lineStyle', '-',  'marker', 'v'), ...
        struct('color', [0.30, 0.75, 0.93], 'lineStyle', '-', 'marker', 'p') ...
    };
end

function ber = process(methodCfg, U, Block_Num, N, C, L, Equal, SNR)
    % 一次完整的仿真链路：发射 -> 信道 -> 接收 -> 计算BER
    [Bits,Symbols0] = Transmitter(U, Block_Num, N, C, methodCfg);
    [H0,Symbols1]   = Channel(Symbols0, L, C, N, Block_Num, SNR);
    Bitsre          = Receiver(U, Block_Num, N, C, Equal, Symbols1, H0, SNR, methodCfg);
    total_bits      = Block_Num * N * log2(U);
    err_bits        = sum(Bits ~= Bitsre);
    ber             = err_bits / total_bits;
end