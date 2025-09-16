clc; clear; close all;

%% —— 基本物理与系统参数 ——    
D2 = 0.5;            % 观测孔径直径 (m)  
wvl = 1e-6;           % 波长 (m)  
k = 2*pi/wvl;       % 波数  
Dz = 50e3;           % 总传播距离 (m)  

%% —— 初始网格与采样设置 ——   
N = 512;                          % 顶点数量  
d1 = 10e-3;                        % 源面采样间距 (m)  
radius = (N*d1)/sqrt(pi);          % 圆盘半径（可改为 D2/2 以匹配口径） 
[pts, TR] = generate_Fibonacci_mesh(N, radius);

% 可视化网格  
figure(1);
triplot(TR);
axis equal;
title('初始 Fibonacci Delaunay 网格');
xlabel('x (m)'); ylabel('y (m)');

%% —— 普通高斯光束参数 ——   
z0 = 0;        % 初始面放在 z=0   
w0 = 0.02;     % 束腰 2 cm

% 计算初始光场（高斯光束）  
U0 = Gaussian_Beam(pts, w0, z0, k);
I0 = abs(U0).^2;
figure(2);
trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I0, 'EdgeColor','none');
view(2); shading interp; axis equal; colorbar;
title('初始高斯光束辐照度分布');

%% 构建传播算子
% 生成初始图拉普拉斯矩阵
L_C = cotangent_Graph_Laplacian(pts,TR);

% 特征分解
opts_eigs = struct('issym',true,'isreal',false);
[Phi, Lambda] = eigs(L_C, N, 'smallestabs', opts_eigs);

%% —— 湍流相关量 —— 
alpha = 5/3;       % Kolmogorov 指数
Cn2   = 1e-16;         % 湍流结构常数 
L0 = 100;              % 外尺度 
l0 = 0.01;             % 内尺度 

r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);   % 平面波 Fried 参量 
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);       % 球面波 Fried 参量 
p = linspace(0,Dz,1e3); 
rytov = 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6)*(p(2)-p(1))); 
    
%% —— 优化求解各相位屏 r0scrn —— 
scr_count = 11;       % 相位屏数量 
n = scr_count; 
A = zeros(2,n); 
alpha_vec = (0:n-1)/(n-1); 
A(1,:) = alpha_vec.^(5/3); 
A(2,:) = (1-alpha_vec).^(5/6).*alpha_vec.^(5/6); 
b = [r0sw^(-5/3); rytov/1.33*(k/Dz)^(5/6)]; 
x0 = (n/3*r0sw*ones(n,1)).^(-5/3); 
x1 = zeros(n,1); 
rmax = 0.1; 
x2 = rmax/1.33*(k/Dz)^(5/6)./A(2,:); 
x2(A(2,:)==0) = 50^(-5/3); 
opts = optimoptions('fmincon','Display','none'); 
[X,~,~,~] = fmincon(@(X) sum((A*X - b).^2), x0, [],[],[],[], x1, x2, [], opts); 
r0scrn = X.^(-3/5); 
r0scrn(isinf(r0scrn)) = 1e6; 

%% 多层分裂步进传播
M_turb = 100;  % KL 模式数（或者设置为你需要的数量）
dz  = Dz / scr_count;  % 每段距离
nreals  = 40;  % 蒙特卡洛次数

% 预分配
U_layers = cell(nreals, scr_count+1); 
epsilon = zeros(nreals, 1);

for real = 1:nreals 
    U_current = U0;  % 初始光场
    U_layers{real, 1} = U0;  
    remeshed_flag = false;  % 标记网格重构

    for layer = 1:scr_count 
        % 半步衍射
        U_current = propagate_HalfStep(U_current, Phi, Lambda, dz/2, k);

        % 生成 KL 相位屏 
        if layer == 1 || remeshed_flag
            rc_step = r0scrn(layer);  % 当前层的 Fried 半径

            % 根据结构函数构造协方差矩阵 L_G
            L_G = construct_Covariance_Laplacian(pts, r0scrn(layer), alpha, 'nugget', 1e-12);

            % 使用协方差拉普拉斯矩阵 L_G 进行特征值分解
            opts_eigs = struct('issym', true, 'isreal', false);
            m_modes = min(400, size(L_G, 1));  % 限制模式数（例如 400）

            try
                [PsiG, DG] = eigs(L_G, m_modes, 'SM', opts_eigs);  % 直接对 L_G 进行特征值分解
            catch
                [Vfull, Dfull] = eig(full(L_G));  % 备用解法
                [dvals, idx] = sort(real(diag(Dfull)), 'ascend');
                idxs = idx(1:m_modes); 
                PsiG = Vfull(:, idxs); 
                DG = diag(dvals(1:m_modes)); 
            end

            % 提取特征值并修正负值
            eigvalsG = real(diag(DG));  
            eigvalsG(eigvalsG < 0) = 0;  % 去掉负特征值
            sigma_n = sqrt(eigvalsG);    % 标准差（KL模式的方差）

            remeshed_flag = false;  % 重网格标志重置
        end

        % 生成实数值的随机高斯系数
        s_n = sigma_n .* randn(length(sigma_n), 1);
        phi = PsiG * s_n;  % 生成相位

        % 叠加湍流相位
        U_current = U_current .* exp(1i * phi);

        % 后半步衍射
        U_current = propagate_HalfStep(U_current, Phi, Lambda, dz/2, k);

        % 自适应网格（每两层一次）
        if mod(layer, 2) == 0
            [pts, U_current, TR] = adaptive_remeshing_sync(pts, TR, U_current, theta_max, theta_min, K_dB);

            % 重建传播算子
            [K, A_vert] = cotangent_K_and_M(pts, TR); 
            M = spdiags(A_vert, 0, size(K, 1), size(K, 1));

            % 更新 KL 模式
            m_modes = min(size(K, 1), min(400, size(K, 1) - 1));
            [Psi, D] = eigs(K, M, m_modes, 'SM', opts_eigs);
            lambda_vec = real(diag(D));
            remeshed_flag = true;  % 设置网格重构标志
            clear PsiG DG L_G
        end

        % 保存当前层的场
        U_layers{real, layer} = U_current;
    end

    % 误差计算（使用 K, M）
    epsilon(real) = compute_error_KM(pts, TR, U_current);
    U_layers{real, scr_count + 1} = U_current;
end