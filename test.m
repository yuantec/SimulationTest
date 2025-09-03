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

%% —— 湍流相关量 —— 
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