%% 半步衍射函数
% function U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
%     %PROPAGATEHALFSTEP   图谱域半步衍射
%     %  U1 = propagateHalfStep(U0,Phi,Lambda,dz,k)
%     %  U0    - 初始场量 (N×1)
%     %  Phi   - 特征向量 (N×M)
%     %  Lambda- 对角矩阵 (M×M)
%     %  dz    - 步长
%     %  k     - 波数
% 
%     % 1) 投影到谱域
%     U_hat = Phi' * U0;
%     % 2) 乘以相位因子
%     phase = exp(1i * diag(Lambda) * (dz/(2*k)));
%     U_hat = U_hat .* phase;
%     % 3) 反变换回空间域
%     U1 = Phi * U_hat;
% end


% function U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
% % PROPAGATE_HALFSTEP  稳健的谱域半步衍射
% %  U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
% %  - U0: N×1 列向量（场）
% %  - Phi: N×M 特征向量（每列一个模式）
% %  - Lambda: M×M 对角矩阵 或 M×1 向量（特征值）
% %  - dz, k: 标量
% %
% %  返回 U1 (N×1)
% 
%     % 确保 U0 是列向量
%     U0 = U0(:);
% 
%     % 投影到谱域
%     U_hat = Phi' * U0;   % 应为 M×1
% 
%     % 处理 Lambda：统一为向量 lam (M×1)
%     if isvector(Lambda)
%         lam = Lambda(:);
%     else
%         % 矩阵时取其对角
%         lam = diag(Lambda);
%     end
% 
%     % 计算按分量的相位因子（确保 lam 是列向量）
%     phase = exp(-1i * lam * (dz/(2*k)));   % M×1
% 
%     % 按元素相乘（M×1 .* M×1 -> M×1）
%     U_hat = U_hat .* phase;
% 
%     % 反变换回空间域（输出 N×1）
%     U1 = Phi * U_hat;
% end

function U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k, M)
% PROPAGATE_HALFSTEP  谱域半步衍射（支持 M 加权投影）
% U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
% U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k, M)
%
% - U0: N×1 复列向量 (field)
% - Phi: N×m 模态矩阵（每列一个模式）
% - Lambda: m×m 对角矩阵 或 m×1 向量（特征值）
% - dz, k: 标量
% - M: 可选，N×N 质量矩阵或长度 N 的顶点权重向量
%
% 说明:
% - 若传入 M 则使用 M 加权内积投影 a = Phi' * (M * U0)；
%   若 Phi 在 M 内积下未正交，函数会对 Phi 做数值稳定的 M-正交化。
% - 若不传 M，函数退回到原先的无权投影行为（向后兼容）。
%
% 返回:
% - U1: N×1 复列向量，传播后场

    % --- 输入校验与标准化 ---
    U0 = U0(:);                 % 确保列向量
    [N, m] = size(Phi);        % N rows, m modes

    % 统一 Lambda 为向量 lam (m x 1)
    if isvector(Lambda)
        lam = Lambda(:);
    else
        lam = diag(Lambda);
    end
    if length(lam) ~= m
        error('长度不匹配: length(Lambda) = %d, but Phi has %d columns.', length(lam), m);
    end

    % 处理可选 M
    useM = (nargin >= 6) && ~isempty(M);
    if useM
        % 若 M 以向量形式给出（顶点质量），变为稀疏对角矩阵
        if isvector(M)
            if length(M) ~= N
                error('如果 M 为向量，其长度必须等于 Phi 的行数 N = %d (现在为 %d).', N, length(M));
            end
            M = sparse(1:N, 1:N, M(:), N, N);
        else
            % 确认矩阵尺寸
            [m1,m2] = size(M);
            if m1~=N || m2~=N
                error('质量矩阵 M 的尺寸必须为 %d x %d，但得到 %d x %d.', N, N, m1, m2);
            end
        end
    end

    % --- 投影到谱域（加权或不加权） ---
    if useM
        % 先计算 Gram 矩阵 B = Phi' * M * Phi，检查是否接近单位阵
        B = Phi' * (M * Phi);
        Bsym = (B + B') / 2;   % 对称化以去除数值噪声
        % 若 B 非常接近 I，则直接使用；否则做 M-正交化
        tol_rel = 1e-8;
        if norm(Bsym - eye(m), 'fro') > tol_rel * max(1, norm(Bsym,'fro'))
            % 做特征分解并构造正交化变换 T
            [V, S] = eig(Bsym);
            s = real(diag(S));
            % 防止除零或负小值：阈值化
            s(s <= eps) = eps;
            T = V * diag(1 ./ sqrt(s)) * V';
            Phi = Phi * T;    % 更新 Phi，使得 Phi' * M * Phi ≈ I
            % 重新计算 Bsym（debug/验证用途，可注释掉）
            % Bsym = (Phi' * (M * Phi) + (Phi' * (M * Phi))')/2;
        end
        % 加权投影
        a = Phi' * (M * U0);   % m x 1
    else
        % 无权投影（向后兼容）
        a = Phi' * U0;
    end

    % --- 传播相位（按模式逐分量相乘） ---
    phase = exp(-1i * lam * (dz/(2*k)));  % m x 1
    a = a .* phase;

    % --- 逆变换回空间域 ---
    U1 = Phi * a;
end
