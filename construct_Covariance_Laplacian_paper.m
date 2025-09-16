function L_G = construct_Covariance_Laplacian_paper(pts, rc, alpha)
% 构造协方差图的拉普拉斯
% Strict paper literal: W(i,j) = (|ri-rj|/rc)^alpha  (Γs as written in paper)
% pts: N×2 coordinates
% rc, alpha: scalars from paper (e.g. alpha=5/3)
% returns sparse Laplacian L_G = D - W (W has zeros on diagonal or keep diag if you prefer)

N = size(pts,1);
if N<2, error('need at least 2 pts'); end

Dmat = squareform(pdist(pts));          % NxN distances
W = (Dmat ./ rc) .^ alpha;              % paper literal Γs

% often graph adjacency sets diag to 0; follow paper choice -> set diagonal to 0
W(1:N+1:end) = 0;

% build Laplacian
deg = sum(W,2);
Dg = spdiags(deg,0,N,N);
L_G = Dg - sparse(W);

% symmetrize & numeric clean
L_G = (L_G + L_G')/2;
rowSum = sum(L_G,2);
if max(abs(rowSum))>1e-10
    L_G = L_G - spdiags(rowSum,0,N,N);
end
end
