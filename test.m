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

%% —— 构建传播算子（刚度 K 与质量 M） *** PATCHED ***
[K, A_vert] = cotangent_K_and_M(pts, TR);   % returns stiffness and vertex areas
M = spdiags(A_vert, 0, size(K,1), size(K,1));

% opts_eigs: don't include issym/isreal to avoid warning when first arg is matrix
opts_eigs = struct();   % *** PATCHED ***

% choose modal truncation (must be < N)
m_modes = min(size(K,1)-1, min(400, size(K,1)-1));
if m_modes < 1
    m_modes = 1;
end

% generalized eigensolve (K psi = lambda M psi) with robust fallback
try
    [Psi, D] = eigs(K, M, m_modes, 'SM', opts_eigs);  % seek smallest magnitude eigenpairs
    lambda_vec = real(diag(D));
catch ME
    warning('generalized eigs(K,M) failed: %s. Falling back to full generalized solve.', E.message);
    % full generalized eigenproblem (dense) - slower
    [Vfull, Dfull] = eig(full(K), full(M));
    dvals = real(diag(Dfull));
    [dvals_sorted, idx] = sort(dvals,'ascend');
    idxs = idx(1:m_modes);
    Psi = Vfull(:, idxs);
    lambda_vec = dvals_sorted(1:m_modes);
end

% M-orthogonality check
orth_err = norm(Psi' * M * Psi - eye(size(Psi,2)));
fprintf('Initial Psi M-orthogonality error = %.3e\n', orth_err);

%% —— 湍流相关物理量 & 分段 Fried 半径 ——   
alpha = 5/3;       % Kolmogorov 指数  
Cn2 = 1e-16;     % 湍流结构常数

r0sw = (0.423*k^2*Cn2*3/8*Dz)^(-3/5);   % 平面波 Fried 参量    
r0pw = (0.423*k^2*Cn2*Dz)^(-3/5);       % 球面波 Fried 参量   
p = linspace(0,Dz,1e3);
rytov = 0.563*k^(7/6)*sum(Cn2*(1-p/Dz).^(5/6).*p.^(5/6)*(p(2)-p(1)));

% 分段优化求解每层 r0scrn
scr_count = 11;               % 相位屏数量    
n = scr_count;
A = zeros(2,n);
alpha_vec = (0:n-1)/(n-1);
A(1,:) = alpha_vec.^(5/3);
A(2,:) = (1-alpha_vec).^(5/6).*alpha_vec.^(5/6);
b = [r0sw^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
x0 = (n/3*r0sw*ones(n,1)).^(-5/3);
x1 = zeros(n,1);
rmax= 0.1;
x2 = rmax/1.33*(k/Dz)^(5/6)./A(2,:);
x2(A(2,:)==0) = 50^(-5/3);

opts_fmin = optimoptions('fmincon','Display','none');
[X,~,~,~] = fmincon(@(X) sum((A*X - b).^2), x0,[],[],[],[], x1, x2,[], opts_fmin);
r0scrn = X.^(-3/5);
r0scrn(isinf(r0scrn)) = 1e6;

%% —— 多层分裂步进传播 ——   
M_turb = 100;                   % KL 模式数  
dz  = Dz/scr_count;             % 每段距离  
nreals  = 40;                   % 蒙特卡洛次数  

% --- 自适应与传播参数 ---  
theta_max = 30;   % degrees  
theta_min = 1;    % degrees  
K_dB = -30;       % dB  

% 预分配  
U_layers = cell(nreals, scr_count+1);
epsilon = zeros(nreals,1);

% 使用 Psi,lambda_vec,M 为传播算子
for real = 1:nreals
    U_current = U0;
    U_layers{real,1} = U0;
    remeshed_flag = false;

    for layer = 1:scr_count
        % ===== 对称分裂步：先半步衍射（使用 M 内积投影） *** PATCHED ***
        U_current = propagate_HalfStep_mass(U_current, Psi, lambda_vec, M, dz/2, k);

        %—— 生成 KL 相位屏 ——   
        if layer==1 || remeshed_flag
            rc_step = r0scrn(layer);
            % 构造 L_G via structure function -> covariance (with small nugget inside)
            L_G = construct_Covariance_Laplacian_structfunc(pts, rc_step, alpha, 'nugget',1e-12);

            % determine requested KL modes (must be < N)
            Ncur = size(K,1);
            nModesReq = min(M_turb, max(1, Ncur-1));

            % add small relative nugget to L_G to improve conditioning
            nugget_rel = 1e-10 * max(1, mean(abs(diag(L_G))));
            L_G_reg = L_G + nugget_rel * speye(size(L_G,1));

            % robust generalized eigendecomposition for L_G_reg
            PsiG = []; DG = [];
            try
                [PsiG, DG] = eigs(L_G_reg, M, nModesReq, 'SM', opts_eigs);
                % validate DG
                if isempty(DG) || ~ismatrix(DG) || size(DG,1) ~= size(DG,2)
                    error('eigs returned unexpected DG (empty or non-square).');
                end
            catch ME
                warning('eigs(L_G_reg,M) failed: %s. Falling back to full generalized solve.', E.message);
                [Vfull, Dfull] = eig(full(L_G_reg), full(M));
                dvals = real(diag(Dfull));
                [dvals_sorted, idx] = sort(dvals,'ascend');
                sel = idx(1:min(nModesReq,length(idx)));
                PsiG = Vfull(:, sel);
                DG = diag(dvals_sorted(1:length(sel)));
            end

            % extract eigenvalues robustly
            if ismatrix(DG) && size(DG,1) == size(DG,2)
                eigvalsG = real(diag(DG));
            elseif isvector(DG)
                eigvalsG = real(DG(:));
            else
                % try conversion
                try
                    eigvalsG = real(diag(full(DG)));
                catch
                    error('Cannot interpret DG as eigenvalue container.');
                end
            end
            eigvalsG(eigvalsG < 0) = 0;
            sigma_n = sqrt(eigvalsG);
            remeshed_flag = false;
        end

        % safety: if sigma_n not defined (very rare), set minimal
        if exist('sigma_n','var') ~= 1 || isempty(sigma_n)
            sigma_n = 1e-12 * ones(max(1, min(10, M_turb)),1);
        end

        % 生成实系数 KL：s_n ~ N(0, sigma_n^2) 并构造相位 phi
        s_n = sigma_n .* randn(length(sigma_n),1);
        phi = PsiG * s_n;
        fprintf('real %d layer %d: std(phi)=%.3e\n', real, layer, std(phi));

        %—— 叠加湍流相位 & 后半步衍射 ——  
        U_current = U_current .* exp(1i * phi);
        U_current = propagate_HalfStep_mass(U_current, Psi, lambda_vec, M, dz/2, k);

        %—— 自适应网格（每2层一次） ——   
        if mod(layer,2)==0
            % 重构网格并同步场 U
            [pts, U_current, TR] = adaptive_remeshing_sync(pts, TR, U_current, theta_max, theta_min, K_dB);

            % 重建 K, M, 和传播谱 Psi, lambda_vec
            [K, A_vert] = cotangent_K_and_M(pts, TR);
            M = spdiags(A_vert, 0, size(K,1), size(K,1));

            % choose new modal truncation (must be < new N)
            newN = size(K,1);
            m_modes = min(newN-1, min(400, newN-1));
            if m_modes < 1, m_modes = 1; end

            try
                [Psi, D] = eigs(K, M, m_modes, 'SM', opts_eigs);
                lambda_vec = real(diag(D));
            catch ME
                warning('eigs(K,M) after remesh failed: %s. Falling back to full generalized solve.', E.message);
                [Vfull, Dfull] = eig(full(K), full(M));
                dvals = real(diag(Dfull));
                [dvals_sorted, idx] = sort(dvals,'ascend');
                idxs = idx(1:m_modes);
                Psi = Vfull(:, idxs);
                lambda_vec = dvals_sorted(1:m_modes);
            end

            remeshed_flag = true;
            clear PsiG DG L_G
        end

        U_layers{real,layer} = U_current;
    end

    % 误差计算（使用 K, M）
    epsilon(real) = compute_error_KM(pts, TR, U_current);
    U_layers{real, scr_count+1} = U_current;
end

%% 可视化示例：第1次试验的初始、中间和最终强度  
for layer = [1, 6, scr_count+1]
    I = abs(U_layers{1,layer}).^2;
    figure;
    trisurf(TR.ConnectivityList, pts(:,1), pts(:,2), I, 'EdgeColor','none');
    view(2); shading interp; axis equal; colorbar;
    if layer==1
        title('试验1：初始场');
    elseif layer==scr_count+1
        title('试验1：最终场');
    else
        title(sprintf('试验1：第 %d 层后', layer-1));
    end
end
