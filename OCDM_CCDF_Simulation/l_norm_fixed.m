function [z_opt, papr_iter] = l_norm_fixed(z, comm_idx, radar_idx, l, iteration)
% z输入波形（频域码元，N×1复数向量）
% comm_idx 通信子载波索引
% radar_idx 雷达子载波索引（优化变量）
% l 范数级次（推荐值10~100，越大越逼近无穷范数）
% iteration 迭代次数
% z_opt 输出优化后的频域波形（N×1复数向量）
% papr_iter 每次迭代的 PAPR (dB)  === 新增输出 ===

N = length(z);
M = 4 * N;

% === 新增：预存储每轮迭代 PAPR ===
papr_iter = zeros(1, iteration);

% 参数检查
if nargin < 5
    error('l_norm_fixed: 需要5个输入参数 (z, comm_idx, radar_idx, l, iteration)');
end

% 索引合法性检查
if any(comm_idx < 1) || any(comm_idx > N)
    error('comm_idx 超出范围 [1, N=%d]', N);
end
if ~isempty(radar_idx) && (any(radar_idx < 1) || any(radar_idx > N))
    error('radar_idx 超出范围 [1, N=%d]', N);
end
if ~isempty(intersect(comm_idx, radar_idx))
    error('comm_idx 与 radar_idx 存在重叠子载波，请检查索引');
end

% 生成4倍过采样FFT矩阵 A (4N×N)
% 对应论文式(3-70): x[n] = (1/√N) Σ X[k]·e^{j2πkn/4N}
% 这里用1/√(4N)归一化以保持Parseval能量守恒
m = (0:N-1)';     % N×1列向量
n = 0:4*N-1;      % 1×4N行向量
exp_term = exp(1j * 2 * pi * m * n / (4*N));  % N×4N
A = (1/sqrt(N)) * exp_term';  % 4N×N，注意论文用正号exp(+j2πkn)

% 保存原始通信码（保持不变）
z_comm = z(comm_idx);

% 初始化：使用输入z作为初始点
z_opt = z;

% 迭代优化
for it = 1:iteration

    % 1. 计算当前频域码对应的4倍过采样时域信号
    % 对应论文式(3-70): x = A * X
    x = A * z_opt;  % 4N×1

    % === 新增：计算并保存当前迭代 PAPR ===
    papr_val = max(abs(x).^2) / mean(abs(x).^2);
    papr_iter(it) = 10*log10(papr_val);  % 转dB存储

    % 2. 计算时域功率 p[n] = |x[n]|^2
    p = abs(x).^2;  % 4N×1

    % 3. 计算l-范数 t = (Σ p_n^l)^(1/l)
    t = sum(p.^l)^(1/l);

    % 4. 计算辅助参数 alpha_n, beta_n
    % 对应论文式(3-104)附近的MM推导
    p_l = p.^l;                  % p_n^l
    p_l_minus_1 = p.^(l-1);      % p_n^(l-1)
    t_minus_p = t - p;           % t - p_n

    % 处理分母为0的情况（数值稳定性）
    denominator = t_minus_p.^2;
    denominator(denominator < eps) = eps;

    % alpha_n = (t^l - p_n^l - l * p_n^(l-1) * (t - p_n)) / (t - p_n)^2
    alpha = (t.^l - p_l - l .* p_l_minus_1 .* t_minus_p) ./ denominator;

    % beta_n = l * p_n^(l-1) - 2 * alpha_n * p_n
    beta = l * p_l_minus_1 - 2 * alpha .* p;

    % 5. 计算gamma = max(2*alpha_n*p_n + beta_n)
    % 对应论文式(3-88)的Lmax
    gamma = max(2 * alpha .* p + beta);

    % 6. 构造对角矩阵（4N×4N）
    D_alpha = diag(alpha);
    D_p = diag(p);
    D_beta = diag(beta);
    I_mat = eye(M);

    % 7. 计算更新方向 Y
    % 对应论文式(3-104): Y = A^H [2*diag(alpha)*diag(p) + diag(beta) - gamma*I] * A * X
    term = 2 * D_alpha * D_p + D_beta - gamma * I_mat;
    Y = - A' * (term * x);  % N×1

    % 8. 闭式解：只取相位（恒模约束）
    % 对应论文式(3-105): X[k] = e^{j*arg(Y[k])}, k ∈ Io
    % 同时保持原始通信码不变
    Y_radar = Y(radar_idx);
    z_opt(radar_idx) = exp(1j * angle(Y_radar));

    % 9. 保持通信码不变（约束条件）
    z_opt(comm_idx) = z_comm;

end

% === 迭代结束后绘制 PAPR 收敛曲线 ===
% figure;
%plot(1:iteration, papr_iter, 'b-o', 'LineWidth',1.5,'MarkerSize',4);
%grid on; grid minor;
%xlabel('迭代次数 Iteration','FontSize',12);
%ylabel('PAPR (dB)','FontSize',12);
%title('l-norm 优化 PAPR 收敛曲线','FontSize',14);

end