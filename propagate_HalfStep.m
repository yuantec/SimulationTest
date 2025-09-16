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


function U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
% PROPAGATE_HALFSTEP  稳健的谱域半步衍射
%  U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
%  - U0: N×1 列向量（场）
%  - Phi: N×M 特征向量（每列一个模式）
%  - Lambda: M×M 对角矩阵 或 M×1 向量（特征值）
%  - dz, k: 标量
%
%  返回 U1 (N×1)

    % 确保 U0 是列向量
    U0 = U0(:);

    % 投影到谱域
    U_hat = Phi' * U0;   % 应为 M×1

    % 处理 Lambda：统一为向量 lam (M×1)
    if isvector(Lambda)
        lam = Lambda(:);
    else
        % 矩阵时取其对角
        lam = diag(Lambda);
    end

    % 计算按分量的相位因子（确保 lam 是列向量）
    phase = exp(-1i * lam * (dz/(2*k)));   % M×1

    % 按元素相乘（M×1 .* M×1 -> M×1）
    U_hat = U_hat .* phase;

    % 反变换回空间域（输出 N×1）
    U1 = Phi * U_hat;
end
