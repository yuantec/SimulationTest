% 基于高斯光的湍流光传输流程验证

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
[pts, TR] = generate_Fibonacci_mesh(N, radius);

%% 初始高斯束
z0 = 0;
w0 = 0.02;
U = Gaussian_Beam(pts, w0, z0, k);

figure(1);
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), abs(U).^2, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar; title('初始强度');

%% 传播算子
L_C = cotangent_Graph_Laplacian(pts, TR);
M_phi = min(300, size(L_C,1)-1);
[Phi, Dlam] = eigs(L_C, M_phi, 'SM');   % 取低频谱模态用于衍射
lambda_vec = real(diag(Dlam));          % M_phi x 1

%% 无湍流传输
U0 = U;
U1 = propagate_HalfStep(U0, Phi, lambda_vec, Dz/2, k);
U2 = propagate_HalfStep(U1, Phi, lambda_vec, Dz/2, k);
figure;
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), abs(U2).^2, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar; title('真空传播光强');

%% 高斯光束在真空中的传播验证
