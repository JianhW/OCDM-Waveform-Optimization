function C_opt = optimize_waveform_freq(c, C_idx, R_idx, varargin)
% OPTIMIZE_WAVEFORM_FREQ  通信-雷达一体化 OFDM 波形优化（输出频域）
%
% 功能：
%   对输入频域码元 c 进行迭代 ICF（削峰+CVX凸优化）优化，
%   在通信子载波 EVM 最小、雷达子载波清零的约束下将 PAPR 压缩到目标阈值，
%   最终返回优化后的 **频域** 码元向量 C_opt（长度 N）。
%
% ─── 调用格式 ───────────────────────────────────────────────────────────
%   C_opt = optimize_waveform_freq(c, C_idx, R_idx)
%   C_opt = optimize_waveform_freq(c, C_idx, R_idx, Name, Value, ...)
%
% ─── 必选输入 ────────────────────────────────────────────────────────────
%   c       - 初始频域码元，(N×1) 复数列向量
%               N = 总子载波数
%   C_idx   - 通信子载波索引，整数向量，范围 1…N
%               EVM 约束仅作用于这些位置
%   R_idx   - 雷达子载波索引，整数向量，范围 1…N
%               优化后该频点幅度置零（可为空 []）
%
% ─── 可选 Name-Value 参数 ─────────────────────────────────────────────────
%   'L'       - 过采样因子（正整数），默认 4
%   'Iter'    - 迭代次数，默认 2
%   'CR'      - 削峰比 = sqrt(10^(PAPR_target_dB/10))
%               默认 sqrt(10^(5/10)) ≈ 1.778（目标 PAPR ≈ 5 dB）
%   'UseCVX'  - 是否使用 CVX 凸优化（true/false），默认 true
%   'A'       - 预计算 IFFT 映射矩阵（LN×N），若不提供则内部生成
%               A(i,k) = (1/√LN) * exp(j·2π·(i-1)·(k-1) / LN)
%   'Verbose' - 是否打印进度（true/false），默认 false
%
% ─── 输出 ─────────────────────────────────────────────────────────────────
%   C_opt   - 优化后的频域码元，(N×1) 复数列向量
%               C_opt(C_idx)：已优化的通信子载波
%               C_opt(R_idx)：强制置零（雷达子载波）
%
% ─── 算法流程 ─────────────────────────────────────────────────────────────
%   初始化：x = ifft([c; 0…0], LN) * √LN
%   For it = 1 : Iter
%     Step 1: T_clip = CR · ‖x‖/√LN；  x_hat = clip(x, T_clip)
%     Step 2: c_in   = fft(x_hat)/√LN 的前 N 点（带内分量）
%     Step 3: T_peak = CR · ‖x_hat‖/√LN
%     Step 4: CVX 求最优滤波器 H（N×1）：
%               min  t
%               s.t. ‖c(C_idx) − H(C_idx)·c_in(C_idx)‖ ≤ ‖c(C_idx)‖·t
%                    H(R_idx)·c_in(R_idx) = 0
%                    |A·(H⊙c_in)| ≤ T_peak
%     Step 5: C_new(1:N) = H⊙c_in；x = ifft(C_new)*√LN
%     Step 6: x = clip(x, T_peak)
%   End
%   输出：C_opt = fft(x)/√LN 的前 N 点；R_idx 处置零
%
% ─── 与 optimize_waveform_cvx.m 的区别 ───────────────────────────────────
%   optimize_waveform_cvx  → 输出时域波形 x_opt  (LN×1)
%   optimize_waveform_freq → 输出频域码元 C_opt  (N×1)   ← 本函数
%
% ─── 示例 ─────────────────────────────────────────────────────────────────
%   N = 128;  L = 4;
%   data = randi([0 3], N, 1);
%   c = pskmod(data, 4, pi/4, 'gray');
%   % ISAC：前 96 个通信子载波，后 32 个雷达子载波
%   C_opt = optimize_waveform_freq(c, 1:96, 97:128, 'L', L, 'Verbose', true);
%
% 版本：v1.0  (2026-03-30)

%% ──────────────── 参数检验与解析 ────────────────────────────────────────
if nargin < 3
    error('optimize_waveform_freq: 至少需要 3 个输入参数 (c, C_idx, R_idx)');
end

% 强制列向量
c     = c(:);
C_idx = C_idx(:)';   % 行向量
R_idx = R_idx(:)';   % 行向量

N = length(c);

% 索引合法性检查
if any(C_idx < 1) || any(C_idx > N)
    error('C_idx 超出范围 [1, N=%d]', N);
end
if ~isempty(R_idx) && (any(R_idx < 1) || any(R_idx > N))
    error('R_idx 超出范围 [1, N=%d]', N);
end
if ~isempty(intersect(C_idx, R_idx))
    error('C_idx 与 R_idx 存在重叠子载波，请检查索引');
end

% 可选参数
p = inputParser;
addParameter(p, 'L',       4,                    @(x) isscalar(x) && x>=1 && mod(x,1)==0);
addParameter(p, 'Iter',    2,                    @(x) isscalar(x) && x>=1 && mod(x,1)==0);
addParameter(p, 'CR',      sqrt(10^(5/10)),      @(x) isscalar(x) && x>0);
addParameter(p, 'UseCVX',  true,                 @islogical);
addParameter(p, 'A',       [],                   @(x) isnumeric(x) || isempty(x));
addParameter(p, 'Verbose', false,                @islogical);
parse(p, varargin{:});

L       = p.Results.L;
iters   = p.Results.Iter;
CR      = p.Results.CR;
use_cvx = p.Results.UseCVX;
A       = p.Results.A;
verbose = p.Results.Verbose;

LN = L * N;

%% ──────────────── 预生成 IFFT 变换矩阵 A ────────────────────────────────
% A(i,k) = (1/√LN) * exp(j·2π·(i-1)·(k-1)/LN)，即 A*c_in ≡ ifft([c_in;0…0])*√LN
% 注意：这里用的是 +j（IFFT 约定），与 fft 的 -j 相反
if isempty(A)
    n_vec = (0:LN-1).';
    k_vec = 0:(N-1);
    A = (1/sqrt(LN)) * exp(1j * 2*pi * (n_vec * k_vec) / LN);  % LN×N
end

if size(A,1) ~= LN || size(A,2) ~= N
    error('A 矩阵维度错误：期望 (%d×%d)，实际 (%d×%d)', LN, N, size(A,1), size(A,2));
end

%% ──────────────── CVX 可用性检测 ────────────────────────────────────────
if use_cvx
    % 尝试加载 CVX
    if exist('cvx_begin', 'file') ~= 2
        if exist('D:\cvx', 'dir')
            addpath(genpath('D:\cvx'));
        end
        try
            cvx_setup; %#ok<NOEFF>
        catch
            use_cvx = false;
            warning('optimize_waveform_freq: CVX 初始化失败，回落到近似求解');
        end
    end

    if use_cvx && exist('cvx_begin', 'file') == 2
        cvx_quiet(true);
        cvx_precision high;
        % 循环尝试 MOSEK → SDPT3，任一可用即停止
        solver_ok = false;
        for try_solver = {'mosek', 'sdpt3'}
            try
                cvx_solver(try_solver{1});
                cvx_begin quiet
                    variable z_probe
                    minimize(z_probe)
                    subject to; z_probe >= 0
                cvx_end
                if strcmpi(cvx_status, 'Solved')
                    if verbose
                        fprintf('[CVX] 使用求解器: %s\n', try_solver{1});
                    end
                    solver_ok = true;
                    break;
                end
            catch
                % 当前求解器不可用，继续尝试下一个
            end
        end
        if ~solver_ok
            use_cvx = false;
            warning('optimize_waveform_freq: CVX 求解器均不可用，回落到近似求解');
        end
    else
        use_cvx = false;
    end
end

%% ──────────────── 时域初始化 ────────────────────────────────────────────
% x = ifft([c; 0…0], LN) * √LN  （√LN 归一化使 ‖x‖² = ‖c‖²）
c_padded = [c; zeros(LN - N, 1)];
x = ifft(c_padded) * sqrt(LN);

%% ──────────────── 主迭代 ────────────────────────────────────────────────
stat = struct('cvx_used', 0, 'cvx_solved', 0, 'fallback', 0, 'approx', 0);

for it = 1:iters

    %--- Step 1: 削峰 ---
    T_clip = CR * norm(x) / sqrt(LN);
    x_hat  = local_clip(x, T_clip);

    %--- Step 2: 提取带内频域分量 ---
    % fft(x_hat)/√LN → 与时域归一化对应
    C_hat = fft(x_hat) / sqrt(LN);
    c_in  = C_hat(1:N);    % N×1

    %--- Step 3: 时域峰值约束阈值 ---
    T_peak = CR * norm(x_hat) / sqrt(LN);

    %--- Step 4: 求最优频域滤波器 H ---
    if use_cvx && exist('cvx_begin', 'file') == 2
        try
            stat.cvx_used = stat.cvx_used + 1;

            if isempty(R_idx)
                %--- 纯通信场景：无雷达子载波清零约束 ---
                cvx_begin quiet
                    variable H(N) complex
                    variable t
                    minimize(t)
                    subject to
                        norm(c(C_idx) - H(C_idx) .* c_in(C_idx), 2) ...
                            <= norm(c(C_idx), 2) * t
                        abs(A * (H .* c_in)) <= T_peak
                cvx_end
            else
                %--- ISAC 场景：含雷达子载波清零约束 ---
                % 将 H(R_idx)·c_in(R_idx) = 0 转写为矩阵形式，
                % 避免在 CVX 块内写 if 语句（CVX 解析器不支持）
                nr  = length(R_idx);
                S_r = zeros(N, nr);
                for ki = 1:nr; S_r(R_idx(ki), ki) = 1; end
                D_radar = diag(c_in(R_idx));   % nr×nr

                cvx_begin quiet
                    variable H(N) complex
                    variable t
                    minimize(t)
                    subject to
                        norm(c(C_idx) - H(C_idx) .* c_in(C_idx), 2) ...
                            <= norm(c(C_idx), 2) * t
                        D_radar * (S_r' * H) == 0
                        abs(A * (H .* c_in)) <= T_peak
                cvx_end
            end

            if strcmpi(cvx_status, 'Solved')
                stat.cvx_solved = stat.cvx_solved + 1;
                H_opt = H;
            else
                if verbose
                    fprintf('[CVX WARN] it=%d, status=%s → feas_project\n', it, cvx_status);
                end
                stat.fallback = stat.fallback + 1;
                H_opt = local_feas_project(ones(N,1), c_in, LN, N, T_peak);
            end

        catch ME
            if verbose
                fprintf('[CVX ERROR] it=%d: %s\n', it, ME.message);
            end
            use_cvx = false;
            stat.approx = stat.approx + 1;
            H_opt = local_feas_project(ones(N,1), c_in, LN, N, T_peak);
        end
    else
        stat.approx = stat.approx + 1;
        H_opt = local_feas_project(ones(N,1), c_in, LN, N, T_peak);
    end

    %--- Step 5: 重建时域波形 ---
    C_new      = zeros(LN, 1);
    C_new(1:N) = H_opt .* c_in;
    % 雷达子载波强制清零（双保险，防止数值误差）
    if ~isempty(R_idx)
        C_new(R_idx) = 0;
    end
    x = ifft(C_new) * sqrt(LN);

    %--- Step 6: 二次轻削峰 ---
    x = local_clip(x, T_peak);
end

%% ──────────────── 提取优化后频域码元（输出） ────────────────────────────
% 对最终时域波形做 FFT，取前 N 点作为频域输出
C_full = fft(x) / sqrt(LN);
C_opt  = C_full(1:N);

% 雷达子载波强制置零（清除数值残差）
if ~isempty(R_idx)
    C_opt(R_idx) = 0;
end

if verbose
    x_check = ifft([C_opt; zeros(LN-N, 1)]) * sqrt(LN);
    papr_db  = 10 * log10(max(abs(x_check).^2) / mean(abs(x_check).^2));
    evm_rms  = norm(C_opt(C_idx) - c(C_idx)) / norm(c(C_idx)) * 100;
    fprintf('[optimize_waveform_freq] 完成 | PAPR=%.2f dB | EVM=%.2f%% | ', ...
        papr_db, evm_rms);
    fprintf('CVX %d/%d Solved | fallback %d | approx %d\n', ...
        stat.cvx_solved, stat.cvx_used, stat.fallback, stat.approx);
end

end  % ← 主函数结束

%% ══════════════════ 辅助函数（局部函数） ══════════════════════════════════

function y = local_clip(x, T)
%LOCAL_CLIP  幅度限幅：|x(i)| > T 时等比缩放到 T
    y   = x;
    mag = abs(x);
    idx = mag > T;
    if any(idx)
        y(idx) = T * x(idx) ./ mag(idx);
    end
end

function H = local_feas_project(H_init, c_in, LN, N, T_peak)
%LOCAL_FEAS_PROJECT  可行化投影（无 CVX 时的回退方案）
%   freq → time → clip → time → freq 一次投影，使时域峰值 ≤ T_peak
    spec = zeros(LN, 1);
    spec(1:N) = H_init .* c_in;
    xt  = ifft(spec) * sqrt(LN);          % 时域重建
    xt  = local_clip(xt, T_peak);         % 削峰
    Cxt = fft(xt) / sqrt(LN);            % 回频域，归一化一致
    Cxt(N+1:end) = 0;                    % 只保留带内 N 点
    H   = Cxt(1:N) ./ (c_in + 1e-12);
end
