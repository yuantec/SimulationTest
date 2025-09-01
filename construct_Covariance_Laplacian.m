function L_G = construct_Covariance_Laplacian(pts, rc, alpha)
% CONSTRUCTCOVARIANCELAPLACIAN 构建协方差拉普拉斯矩阵 L_G
%   L_G = constructCovarianceLaplacian(pts, rc, alpha)
%   输入:
%       pts   - N×2 网格顶点坐标 [x, y]
%       rc    - 波前相干半径
%       alpha - 相位协方差函数指数 (论文中 Γs(r1,r2) = (|r1-r2|/rc)^α)
%   输出:
%       L_G   - N×N 协方差拉普拉斯矩阵

    % 1) 顶点数 N
    N = size(pts,1);
    % 2) 计算顶点对两两欧氏距离
    Dmat = squareform(pdist(pts));

    % 3) 构造协方差权重矩阵 WΓ
    %    Γs(r_i,r_j) = (|r_i-r_j| / rc)^alpha
    W = exp( - (Dmat ./ rc) .^ alpha );   % r 越大，W 越小
    scale = trace(W)/N;
    if scale > 0
        W = W ./ scale;
    end
    % 对角元素置零（自相关项不参与）
    W(1:N+1:end) = 0;

    % 4) 度矩阵 Dg
    Dg = diag(sum(W,2));

    % 5) 协方差拉普拉斯 L_G = Dg - W
    L_G = Dg - W;

    % 6) 验证（对称、行和为零、半正定）
    assert(issymmetric(L_G), 'L_G 非对称');
    assert(all(abs(sum(L_G,2))<1e-10), 'L_G 行和非零');
    % 半正定性可在外部特征分解时验证
end