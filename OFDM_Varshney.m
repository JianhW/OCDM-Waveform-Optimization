function z_opt = OFDM_Varshney(z, idx_comm, idx_radar, varargin)
% OPTIMIZE_OFDM_PAPR  基于论文 Algorithm 1 的 OFDM ISAC PAPR 优化
%
% 功能说明：
%   固定通信子载波 z(idx_comm) 不变，仅优化雷达子载波 z(idx_radar) 的相位，
%   以最小化时域 PAPR。算法在每次迭代中构造 SOCP 问题（CVX + MOSEK），
%   若求解失败则回退到 one-hot 最优单点策略，最后对雷达子载波取负相位更新。
%
% ─────────────────────────────────────────────────────────────────────────────
% 输入参数
%   z         (N×1 复向量, 必填)
%             初始全频域 OFDM 码（通信 + 雷达子载波均已填入）。
%             要求所有子载波满足单位模 |z[k]| = 1（unimodular）。
%             其中通信子载波 z(idx_comm) 在优化过程中保持不变。
%
%   idx_comm  (Nc×1 整数向量, 必填)
%             通信子载波在 z 中的索引（1-based），即固定子载波的位置。
%
%   idx_radar (Nr×1 整数向量, 必填)
%             雷达子载波在 z 中的索引（1-based），即待优化子载波的位置。
%
% 可选名值对参数（Name-Value）
%   'maxIt'     最大迭代次数（默认 10）
%   'epsStop'   PAPR 相对变化收敛阈值（默认 1e-5）
%   'tol_mag'   相位保护幅度下限，|q|<tol 时继承旧相位（默认 1e-12）
%   'lambda2'   L2 微正则化系数，稳定数值（默认 1e-10）
%   'verbose'   是否打印每轮迭代信息（默认 false）
%
% 输出参数
%   z_opt     (N×1 复向量)
%             优化后的全频域 OFDM 码（通信子载波保持不变，雷达子载波已优化）。
%             满足 |z_opt[k]| = 1（所有子载波单位模）。
%
% ─────────────────────────────────────────────────────────────────────────────
% 使用示例
%
%   % 基本用法（使用默认参数）
%   N = 128; Nc = round(0.3*N); Nr = N - Nc;
%   rand_idx = randperm(N);
%   idx_comm  = sort(rand_idx(1:Nc)).';
%   idx_radar = sort(rand_idx(Nc+1:end)).';
%   zc = exp(1j * pi/4 * randi([0 3], Nc, 1));  % QPSK，单位模
%   zr = exp(1j * 2*pi * rand(Nr, 1));           % 雷达子载波，随机初始相位
%   z0 = complex(zeros(N,1));
%   z0(idx_comm)  = zc;
%   z0(idx_radar) = zr;
%   z_opt = optimize_ofdm_papr(z0, idx_comm, idx_radar);
%
%   % 自定义参数
%   z_opt = optimize_ofdm_papr(z0, idx_comm, idx_radar, ...
%               'maxIt', 20, 'epsStop', 1e-6, 'verbose', true);
%
% ─────────────────────────────────────────────────────────────────────────────
% 算法流程（对应论文 Algorithm 1）
%   0. 预处理：强制 z 单位模，提取基本参数
%   for t = 1 : maxIt
%     1. 构造 Q 矩阵和 c 向量（论文式 13）
%        s  = IFFT(z) * sqrt(N)
%        Q  = 2 * (Fmat .* repmat(s.',N,1) - repmat(z,1,N))
%        c  = real(N - 0.5 * (z^H * Q))
%     2. CVX 求解 SOCP（论文式 21）
%        min   ||Qr*w||_1 - c*w + lambda2*||w||^2
%        s.t.  w >= 0, sum(w) = 1
%        （失败时回退：w = one-hot at argmax(Re(z^H Q) + c)）
%     3. 相位更新（论文式 22）
%        q_w = Q * w
%        zr_new = -exp(j * angle(q_w(idx_radar)))
%        （|q_w| < tol 时继承旧相位）
%     4. 重建 z，强制单位模
%     5. 收敛判断（PAPR 相对变化 < epsStop 则提前退出）
%   end
%
% ─────────────────────────────────────────────────────────────────────────────
% 依赖项
%   - MATLAB Signal Processing Toolbox（ifft, dftmtx 等）
%   - CVX（http://cvxr.com/cvx/）及 MOSEK 求解器（推荐）
%     若 CVX/MOSEK 不可用，算法将自动使用回退策略（one-hot），不报错。
%
% 版本：2026-03-29
% ─────────────────────────────────────────────────────────────────────────────

%% ========== 0. 解析输入参数 ==========
p = inputParser;
addRequired(p, 'z',         @(x) isvector(x) && ~isreal(x) || isnumeric(x));
addRequired(p, 'idx_comm',  @(x) isvector(x) && isnumeric(x));
addRequired(p, 'idx_radar', @(x) isvector(x) && isnumeric(x));
addParameter(p, 'maxIt',    10,    @(x) isscalar(x) && x > 0);
addParameter(p, 'epsStop',  1e-5,  @(x) isscalar(x) && x > 0);
addParameter(p, 'tol_mag',  1e-12, @(x) isscalar(x) && x >= 0);
addParameter(p, 'lambda2',  1e-10, @(x) isscalar(x) && x >= 0);
addParameter(p, 'verbose',  false, @islogical);
parse(p, z, idx_comm, idx_radar, varargin{:});

maxIt    = p.Results.maxIt;
epsStop  = p.Results.epsStop;
tol_mag  = p.Results.tol_mag;
lambda2  = p.Results.lambda2;
verbose  = p.Results.verbose;

%% ========== 1. 预处理 ==========
z         = z(:);           % 确保列向量
idx_comm  = idx_comm(:);
idx_radar = idx_radar(:);
N         = numel(z);
Nr        = numel(idx_radar);

% 提取固定的通信子载波（不因迭代改变）
zc = z(idx_comm);

% 强制初始 z 单位模（避免外部传入非单位模符号导致数值问题）
valid_mask = abs(z) >= tol_mag;
z(valid_mask)  = z(valid_mask) ./ abs(z(valid_mask));
z(~valid_mask) = exp(1j * 2*pi * rand(sum(~valid_mask), 1));  % 极小幅度随机重置

% 单位化 DFT 算子：x = F^H z = ifft(z)*sqrt(N)
FH   = @(v) ifft(v(:)) * sqrt(N);
Fmat = dftmtx(N) / sqrt(N);    % N×N DFT 矩阵（仅用于构造 Q）

% PAPR 工具
papr_fun = @(x) max(abs(x).^2) / mean(abs(x).^2);

%% ========== 2. 迭代优化 ==========
x         = FH(z);
PAPR_prev = papr_fun(x);

for t = 1 : maxIt
    % ----- (1) 构造 Q 矩阵和 c 向量 -----
    s  = FH(z);                                   % N×1 时域向量
    FS = Fmat .* repmat(s.', N, 1);               % N×N，第(n,k)列 = F_{nk} * s_n
    Q  = 2 * (FS - repmat(z, 1, N));              % N×N，第 k 列 = q_k（梯度上界方向）
    c  = real(N - 0.5 * (conj(z).' * Q));         % 1×N，目标函数上界截距
    Qr = Q(idx_radar, :);                         % Nr×N，仅保留雷达子载波行

    % ----- (2) CVX 求解 SOCP（论文式 21） -----
    %   min_{w}   ||Qr * w||_1 - c * w + lambda2 * ||w||_2^2
    %   s.t.      w >= 0, sum(w) = 1
    ok = false;
    w  = [];
    try
        cvx_begin quiet
            cvx_precision best
            variable w(N) nonnegative
            minimize( norm(Qr * w, 1) - c * w + lambda2 * sum_square(w) )
            subject to
                sum(w) == 1
        cvx_end
        ok = (strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')) ...
             && all(isfinite(w));
    catch
        ok = false;
    end

    % ----- 回退策略（CVX 失败时）：one-hot 权重 -----
    if ~ok
        % 选择使目标函数上界最"紧"的单点 n*
        score        = real(conj(z).' * Q) + c;   % 1×N
        [~, nstar]   = max(score);
        w            = zeros(N, 1);
        w(nstar)     = 1;
    end

    % ----- (3) 雷达子载波相位更新（论文式 22） -----
    q_w   = Q * w;                                % N×1，等效梯度方向
    q_w_r = q_w(idx_radar);                       % Nr×1，仅雷达部分

    % 相位保护：|q| 极小时继承旧相位，避免 angle(0) 产生 NaN
    keep_mask = abs(q_w_r) < tol_mag | ~isfinite(q_w_r);
    zr_new    = -exp(1j * angle(q_w_r));          % 论文式(22)，取负相位
    if any(keep_mask)
        zr_new(keep_mask) = z(idx_radar(keep_mask));
    end

    % ----- (4) 重建 z，强制单位模 -----
    z              = complex(zeros(N, 1));
    z(idx_comm)    = zc;                          % 通信子载波保持不变
    z(idx_radar)   = zr_new;                      % 更新雷达子载波
    % 强制所有子载波单位模（防止数值误差累积）
    valid_now      = abs(z) >= tol_mag;
    z(valid_now)   = z(valid_now) ./ abs(z(valid_now));
    z(~valid_now)  = exp(1j * 2*pi * rand(sum(~valid_now), 1));

    % ----- (5) 收敛判断 -----
    x         = FH(z);
    PAPR_curr = papr_fun(x);
    rel       = abs(PAPR_curr - PAPR_prev) / max(PAPR_prev, eps);

    if verbose
        fprintf('  [OFDM PAPR opt] iter %3d/%d: PAPR = %.4f dB, rel_err = %.2e%s\n', ...
            t, maxIt, 10*log10(PAPR_curr), rel, ...
            iif(~ok, ' [fallback]', ''));
    end

    if rel <= epsStop
        break;
    end
    PAPR_prev = PAPR_curr;
end

%% ========== 3. 输出 ==========
z_opt = z;

end  % function optimize_ofdm_papr


%% ========== 局部辅助函数 ==========
function out = iif(cond, a, b)
% 内联三元运算符（避免与用户环境中的 ternary 冲突）
    if cond
        out = a;
    else
        out = b;
    end
end
