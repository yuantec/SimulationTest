% test_splitstep_paperLiteral.m
% 严格按论文文字：W = Gamma_s = (|r|/rc)^alpha -> L_G = D - W
% 单次 split-step 传播（无自适应），简洁实现

clc; clear; close all;

%% 基本物理与网格参数
D2 = 0.5;           % 观测口径 (m)
wvl = 1e-6;         % 波长 (m)
k = 2*pi/wvl;
Dz = 50e3;

%% 网格
N = 512;
d1 = 10e-3;
radius = (N*d1)/sqrt(pi);
[pts, TR] = generate_Fibonacci_mesh(N, radius);  % 你已有该函数

%% 初始高斯束
z0 = 0;
w0 = 0.02;
U = Gaussian_Beam(pts, w0, z0, k);

figure(1);
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), abs(U).^2, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar; title('初始强度');

%% 传播算子：余切图 Laplacian L_C -> 用于谱域衍射
L_C = cotangent_Graph_Laplacian(pts, TR);    % 你已有此函数或用之前版本
M_phi = min(300, size(L_C,1)-1);
[Phi, Dlam] = eigs(L_C, M_phi, 'SM');   % 取低频谱模态用于衍射
lambda_vec = real(diag(Dlam));          % M_phi x 1

%% 湍流 / 相位屏参数（按论文）
alpha = 5/3;
Cn2 = 1e-16;
scr_count = 11;
dz = Dz / scr_count;
M_turb = 100;    % 取相位屏模态数（小于 N）
nreals = 1;

% 分段 r0scrn（沿用你原来计算）
p = linspace(0,Dz,1e3);
r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);
rytov = 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6)*(p(2)-p(1)));
n = scr_count;
alpha_vec = (0:n-1)/(n-1);
A = zeros(2,n);
A(1,:) = alpha_vec.^(5/3);
A(2,:) = (1-alpha_vec).^(5/6).*alpha_vec.^(5/6);
b = [r0sw^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
x0 = (n/3*r0sw*ones(n,1)).^(-5/3);
x1 = zeros(n,1);
rmax = 0.1;
x2 = rmax/1.33*(k/Dz)^(5/6)./A(2,:);
x2(A(2,:)==0) = 50^(-5/3);
opts = optimoptions('fmincon','Display','none');
X = fmincon(@(X) sum((A*X - b).^2), x0, [],[],[],[], x1, x2, [], opts);
r0scrn = X.^(-3/5); r0scrn(isinf(r0scrn)) = 1e6;

%% 逐层 split-step（无自适应）
for layer = 1:scr_count
    % 半步谱域衍射
    U = propagate_HalfStep(U, Phi, lambda_vec, dz/2, k);

    % 按论文字面构造 L_G (W = Gamma_s)
    rc = r0scrn(layer);
    L_G = construct_Covariance_Laplacian_paper(pts, rc, alpha);  % 论文字面

    % 对 L_G 做特征分解并采样模态（论文字面做法）
    Mreq = min(M_turb, size(L_G,1)-1);
    [PsiG, DG] = eigs(L_G, Mreq, 'LM');   % 取“largest magnitude”按论文字面（注意可能有负值）
    eigvals = real(diag(DG));
    % 为保证后续数值稳定，仅取非负分量（你要求简单实现，我不做复杂分支）
    eigvals(eigvals < 0) = 0;
    sigma_modes = sqrt(eigvals);

    % 采样系数并合成相位（实相位）
    coeffs = sigma_modes .* randn(length(sigma_modes),1);
    phi = real(PsiG * coeffs);

    fprintf('layer %2d: std(phi)=%.3e\n', layer, std(phi));

    % 叠加相位屏并后半步衍射
    U = U .* exp(1i * phi);
    U = propagate_HalfStep(U, Phi, lambda_vec, dz/2, k);
end
U11111 = abs(U).^2;
%% 显示最终强度
figure;
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), abs(U).^2, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar; title('最终强度（论文字面 W=Gamma_s）');
