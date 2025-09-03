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

