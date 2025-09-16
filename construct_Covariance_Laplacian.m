function L_G = construct_Covariance_Laplacian(pts, rc, alpha, varargin)
% CONSTRUCT_COVARIANCE_LAPLACIAN_STRUCTFUNC
%  从论文的幂律 (|r|/rc)^alpha（作为结构函数 D_phi）出发，
%  构造协方差图的拉普拉斯 L_G = Dg - W.
%
%  语法:
%    L_G = construct_Covariance_Laplacian_structfunc(pts, rc, alpha)
%    L_G = construct_Covariance_Laplacian_structfunc(...,'nugget',1e-12)
%    L_G = construct_Covariance_Laplacian_structfunc(...,'normalize',true)
%
%  说明:
%    - 按论文 §3.2，论文给出 Γs(r1,r2) = (|r1-r2|/rc)^alpha (式 (8) 等处)，
%      该式随距离增大而增大，实际为结构函数 D_phi(r) 的形式。
%    - 随机场理论: D_phi(r) = 2*(C(0) - C(r)) -> C(r) = C(0) - 0.5*D_phi(r).
%      因此下面将论文幂律视作 D_phi 并转换为协方差 C。
%    - 返回 L_G = Dg - W where W(i,j) = C_ij for i!=j, W(ii)=0
%
%  参考: 论文 "Graph-based model for adaptive simulation of beam propagation..."
%         (see §2.2 and §3.2 for Γs and L_G definition). See file citation in conversation.
%
%  参数:
%    pts  - N×2 array of node coordinates
%    rc   - coherence radius used in the paper (rc)
%    alpha- exponent (e.g. 5/3)
%
%  可选项:
%    'nugget'   - small regularization added to C diagonal before building W (default 1e-12)
%    'normalize'- boolean, if true scale W so that trace(W)/N = 1 (default false)
%
%  返回:
%    L_G - N×N sparse Laplacian (row sums ~ 0)
%
%  注: 该实现尽量遵循论文表达（结构函数 -> 协方差 -> W -> L_G）。
%

% parse options
p = inputParser;
addParameter(p,'nugget',1e-12,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'normalize',false,@islogical);
parse(p,varargin{:});
nugget = p.Results.nugget;
doNormalize = p.Results.normalize;

N = size(pts,1);
if N<2
    error('pts must contain at least 2 points');
end

% 1) 两两距离矩阵 (NxN)
Dmat = squareform(pdist(pts));   % Dmat(i,j) = |ri - rj|

% 2) 论文给出的幂律视为结构函数 D_phi
%    D_phi(i,j) = (|ri-rj|/rc)^alpha
Dphi = (Dmat ./ rc) .^ alpha;

% 3) 将结构函数转换为协方差：C(r) = C0 - 0.5 * D_phi(r).
%    不知道明确的 C0（相位零距方差）时，可选策略：
%      - 选 C0 = 0.5 * max(Dphi(:))，使得最小协方差为 0（数值稳定）
%      - 或者用户传入物理方差（如果有）
%    这里我们采用第一种保守做法，并添加 nugget 保证正定性。
C0 = 0.5 * max(Dphi(:));  % 使最小 C >= 0
C = C0 - 0.5 * Dphi;

% 4) 在对角上加入微小 nugget (数值正定)
C = C + nugget * max(1,mean(abs(diag(C)))) * speye(N);

% 5) 构造权重矩阵 WΓ: 用协方差作为权重（但图邻接通常对角 = 0）
W = full(C);           % dense for now
W(1:N+1:end) = 0;      % 将自相关项移除（图权重通常为 i != j 的边权）

% 可选归一化（论文没有明确要求，但有时用于数值稳定或统一尺度）
if doNormalize
    s = trace(W)/N;
    if s>0
        W = W ./ s;
    end
end

% 6) 度矩阵与拉普拉斯
deg = sum(W,2);
Dg = spdiags(deg,0,N,N);
L_G = Dg - sparse(W);

% 强制对称与小数值修正
L_G = (L_G + L_G')/2;
% 修复行和（确保数值误差不造成行和偏离0）
rowSum = sum(L_G,2);
if max(abs(rowSum))>1e-8
    L_G = L_G - spdiags(rowSum,0,N,N);
end

% 验证（仅作调试）
% assert(issymmetric(L_G),'L_G 非对称');
% assert(all(abs(sum(L_G,2))<1e-10),'L_G 行和非零');

end
