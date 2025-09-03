%% 半步衍射函数
function U1 = propagate_HalfStep(U0, Phi, Lambda, dz, k)
    %PROPAGATEHALFSTEP   图谱域半步衍射
    %  U1 = propagateHalfStep(U0,Phi,Lambda,dz,k)
    %  U0    - 初始场量 (N×1)
    %  Phi   - 特征向量 (N×M)
    %  Lambda- 对角矩阵 (M×M)
    %  dz    - 步长
    %  k     - 波数

    % 1) 投影到谱域
    U_hat = Phi' * U0;
    % 2) 乘以相位因子
    phase = exp(1i * diag(Lambda) * (dz/(2*k)));
    U_hat = U_hat .* phase;
    % 3) 反变换回空间域
    U1 = Phi * U_hat;
end