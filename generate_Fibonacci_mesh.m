%% Delauney三角网格初始生成
function [pts, TR] = generate_Fibonacci_mesh(N, R)
%生成斐波那契螺旋点并进行 Delaunay 三角剖分
%   [pts, TR] = generateFibonacciDelaunay(N, R)
%   输入:
%       N - 顶点数量
%       R - 圆盘半径
%   输出:
%       pts - N×2 数组，包含所有顶点坐标 [x, y]
%       TR  - delaunayTriangulation 对象

    % 黄金比例与黄金角
    %N = 1024;
    %R = 1.0;
    phi = (1 + sqrt(5)) / 2;
    goldenAngle = 2*pi*(1 - 1/phi);

    % 预分配
    pts = zeros(N, 2);
    for k = 0:(N-1)
        % 半径 r_k = R * sqrt((k+0.5)/N)
        r = R * sqrt((k + 0.5) / N);
        % 角度 theta_k = k * goldenAngle
        theta = k * goldenAngle;
        % 转为笛卡尔坐标
        pts(k+1,1) = r * cos(theta);
        pts(k+1,2) = r * sin(theta);
    end

    % Delaunay 三角剖分
    TR = delaunayTriangulation(pts(:,1), pts(:,2));

    % 可视化检查（可选）
    % triplot(TR);
    % axis equal;
    % title(sprintf('Fibonacci Delaunay Mesh: N=%d, R=%.2f', N, R));
end
