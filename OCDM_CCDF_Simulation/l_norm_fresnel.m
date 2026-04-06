function [z_opt, papr_iter] = l_norm_fresnel(z, comm_idx, radar_idx, l, iteration, varargin)
% l_norm_fresnel: 基于Fresnel变换的LNCA算法
% 将原始FFT基替换为Fresnel-chirp基
%
% 输入参数:
%   z - 输入波形（频域码元，N×1复数向量）
%   comm_idx - 通信子载波索引
%   radar_idx - 雷达子载波索引（优化变量）
%   l - 范数级次（推荐值10~100，越大越逼近无穷范数）
%   iteration - 迭代次数
%   varargin{1} - 可选：Fresnel参数 [lambda, delta_x, z_distance]
%               lambda: 波长 (默认 0.03m，对应30GHz)
%               delta_x: 采样间隔 (默认 0.5m)
%               z_distance: 传播距离 (默认 100m)
%
% 输出参数:
%   z_opt - 输出优化后的频域波形（N×1复数向量）
%   papr_iter - 每次迭代的 PAPR (dB)
%
% Fresnel变换核 (含标准常数相位因子):
%   A(n,k) = exp(-jπ/4) / √N · exp(j*chirp_rate*(n-k)^2)
%   其中 chirp_rate = π·Δx² / (λ·z)
%
%   注: exp(-jπ/4) = (1-j)/√2，来自Fresnel积分渐近展开的常数相位

N = length(z);
M = 4 * N;  % 4倍过采样

% === 新增：预存储每轮迭代 PAPR ===
papr_iter = zeros(1, iteration);

% 参数检查
if nargin < 5
    error('l_norm_fresnel: 需要5个输入参数 (z, comm_idx, radar_idx, l, iteration)');
end

% 解析Fresnel参数
if nargin >= 6 && ~isempty(varargin{1})
    lambda = varargin{1}(1);
    delta_x = varargin{1}(2);
    z_dist = varargin{1}(3);
else
    % 默认参数（30GHz雷达场景）
    lambda = 0.03;    % 波长 0.03m (对应 f=10GHz)
    delta_x = 0.5;    % 采样间隔 0.5m
    z_dist = 100;     % 传播距离 100m
end

% 计算Fresnel chirp系数
% chirp_rate = pi / (lambda * z) * (delta_x^2)
% 对应Fresnel变换核: exp(j * chirp_rate * n^2)
chirp_rate = pi * (delta_x^2) / (lambda * z_dist);

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

% === 生成Fresnel变换矩阵 A_fresnel (4N×N) ===
% 标准Fresnel变换核:
%   A_fresnel(n,k) = C · exp(j*chirp_rate*(n - k)^2)
%
% 其中常数相位因子 C = exp(-jπ/4) = (1-j)/√2
% 来源: ∫exp(jπx²)dx ∝ exp(jπ/4), 对应Fresnel数渐近展开
%
% 对应论文式(3-70)的FFT核: exp(j*2π*k*n/4N)
% Fresnel变换是时不变系统，核函数为chirp信号

m = (0:N-1)';     % N×1列向量 (频域索引 k)
n = 0:4*N-1;      % 1×4N行向量 (时域采样 n)

% Fresnel核: exp(j * chirp_rate * (n - m)^2)
chirp_term = exp(1j * chirp_rate * (n(:).^2 - 2*n(:).*m(:)' + m(:).^2'));

% 标准Fresnel变换常数相位因子: exp(-jπ/4)
% 等价形式: (1 - 1j) / sqrt(2) = exp(-jπ/4)
fresnel_const = exp(-1j * pi / 4);  % = (1-j)/√2

% 完整Fresnel变换矩阵（含常数相位因子）
A_fresnel = fresnel_const * (1/sqrt(N)) * chirp_term;  % 4N×N

% 注: 也可使用展开形式的chirp-Z变换（数学等价）:
% n_chirp = exp(1j * chirp_rate * n(:).^2);       % 4N×1 (输出chirp调制)
% k_chirp = exp(-1j * chirp_rate * m(:).^2);      % N×1 (输入chirp调制)
% cross_chirp = exp(1j * 2 * chirp_rate * n(:) * m(:)');  % 4N×N (交叉chirp)
% A_chirpZ = fresnel_const * (1/sqrt(N)) * (n_chirp * k_chirp.') .* cross_chirp;

% 保存原始通信码（保持不变）
z_comm = z(comm_idx);

% 初始化：使用输入z作为初始点
z_opt = z;

% 迭代优化
for it = 1:iteration

    % 1. 计算当前频域码对应的4倍过采样时域信号（使用Fresnel变换）
    % 对应论文式(3-70)的Fresnel版本: x = A_fresnel * X
    x = A_fresnel * z_opt;  % 4N×1

    % === 计算并保存当前迭代 PAPR ===
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
    % 对应论文式(3-104)的Fresnel版本
    % Y = A_fresnel^H [2*diag(alpha)*diag(p) + diag(beta) - gamma*I] * A_fresnel * X
    term = 2 * D_alpha * D_p + D_beta - gamma * I_mat;
    Y = -A_fresnel' * (term * x);  % N×1

    % 8. 闭式解：只取相位（恒模约束）
    % 对应论文式(3-105): X[k] = e^{j*arg(Y[k])}, k ∈ Io
    % 同时保持原始通信码不变
    Y_radar = Y(radar_idx);
    z_opt(radar_idx) = Y_radar;

    % 9. 保持通信码不变（约束条件）
    z_opt(comm_idx) = z_comm;

end

end
