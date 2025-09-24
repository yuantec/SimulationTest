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

%% 测试
L_C = cotangent_Graph_Laplacian(pts, TR);   %刚度矩阵（余切拉普拉斯矩阵）
M = vertex_mass_matrix(pts, TR);  % 质量矩阵
M_phi = min(451, size(L_C,1)-1);
opts.isreal=true;
opts.issym=true;
U0 = U;
% 假设你已经有 L_C, M, U0, opts
% 先算较大的 M_phi_all（尽量接近 N，但注意内存/时间）
M_phi_all = min(511, size(L_C,1)-1);   % N=512 的情形
[Phi_all, Dlam_all] = eigs(L_C, M, M_phi_all, 'SM', opts);
lambda_all = real(diag(Dlam_all));
a_all = Phi_all' * (M * U0);            % 所有模态系数
modal_energy = abs(a_all).^2;
cum_energy = cumsum(modal_energy) / sum(modal_energy);

figure;
plot(1:length(cum_energy), cum_energy, '-o');
xlabel('mode index'); ylabel('cumulative modal energy fraction');
grid on;

% 绘制模态系数幅值（谱）以观察是否有长尾噪声
figure;
semilogy(1:length(modal_energy), modal_energy, '.-');
xlabel('mode index'); ylabel('|a_i|^2 (log scale)');
grid on;

% 可视化若干高阶模式看看它们是否是噪声
for idx = [10, 30, 60, 120, 300, min(480,length(lambda_all))]
    figure;
    trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), real(Phi_all(:,idx)), 'EdgeColor','none');
    view(2); shading interp; axis equal; colorbar; title(sprintf('mode %d, lambda=%.3e', idx, lambda_all(idx)));
end


%% 传播算子
L_C = cotangent_Graph_Laplacian(pts, TR);   %刚度矩阵（余切拉普拉斯矩阵）
M = vertex_mass_matrix(pts, TR);  % 质量矩阵
M_phi = min(100, size(L_C,1)-1);
opts.isreal=true;
opts.issym=true;
[Phi, Dlam] = eigs(L_C, M, M_phi, 'SM', opts);   % L * phi = lambda * M * phi
lambda_vec = real(diag(Dlam));          % M_phi x 1

%% 相关验证
U0 = U;
% 1) mass matrix 是否存在（若没有先算 M）


% 2) 广义正交检验
orthErr = norm(Phi' * (M * Phi) - eye(size(Phi,2)), 'fro');
fprintf('Phi^T * M * Phi - I Frobenius = %.3e\n', orthErr);

% 3) lambda 范围
fprintf('lambda min/max = %.6e / %.6e\n', min(lambda_vec), max(lambda_vec));

% 4) 能量守恒（单步传播前后）
E_before = real(U0' * (M * U0));
Utest = propagate_HalfStep(U0, Phi, lambda_vec, Dz, k, M); % 用你的传播或我给的函数
E_after = real(Utest' * (M * Utest));
fprintf('能量比 E_after / E_before = %.6f\n', E_after / E_before);

% 5) 单一模态传播检验（第 i 个特征向量只应累相位，模态幅值不变）
i = 5;
phi_i = Phi(:,i);
a_i = phi_i' * (M * phi_i); % 应近似 1
phi_i_propag = propagate_HalfStep(phi_i, Phi, lambda_vec, Dz, k, M);
fprintf('mode %d relative energy change = %.3e\n', i, real(phi_i_propag' * (M * phi_i_propag) / a_i - 1));

% 6) 跟解析高斯的 L2 误差（在同一 z）
U_ana = Gaussian_Beam(pts, w0, Dz, k);
L2rel = norm(abs(Utest) - abs(U_ana)) / norm(abs(U_ana));
fprintf('强度 L2 相对误差 = %.3e\n', L2rel);

% 计算投影系数和能量
a = Phi' * (M * U0);              % 模态系数
E_modal = real(a' * a);           % 投影到所选模态的能量（因为 Phi^T M Phi = I）
E_total = real(U0' * (M * U0));   % 原始能量
fprintf('投影模态能量 E_modal = %.6e, 原始能量 E_total = %.6e, 占比 = %.4f\n', E_modal, E_total, E_modal/E_total);

% 计算残差能量（在正交子空间里的能量）
U_rec = Phi * a;
r = U0 - U_rec;
E_res = real(r' * (M * r));
fprintf('残差能量 E_res = %.6e, 占比 = %.4f\n', E_res, E_res/E_total);

% 若残差占比很大，说明需要更多模态或更细网格

%% 无湍流传输
U0 = U;
U1 = propagate_HalfStep(U0, Phi, lambda_vec, Dz/2, k, M);
U2 = propagate_HalfStep(U1, Phi, lambda_vec, Dz/2, k, M);
figure;
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), abs(U2).^2, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar; title('真空传播光强');

%% 高斯光束理论效果
U_ana = Gaussian_Beam(pts, w0, Dz, k);

figure;
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), abs(U_ana).^2, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar; title('真空传播光强');

