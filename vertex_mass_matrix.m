function M = vertex_mass_matrix(pts, TR)
% 计算顶点质量矩阵（lumped mass），每个三角面积平均分配到三个顶点
% 输入:
%   pts - N x 2 顶点坐标
%   TR  - delaunayTriangulation 对象
% 输出:
%   M   - N x N 稀疏对角质量矩阵

T = TR.ConnectivityList;   % Mtri x 3
N = size(pts,1);
vertexArea = zeros(N,1);

% 遍历每个三角形，累加面积到三个顶点（简单可靠）
for t = 1:size(T,1)
    tri = T(t,:);                   % 三个顶点索引
    p1 = pts(tri(1),:);
    p2 = pts(tri(2),:);
    p3 = pts(tri(3),:);
    v1 = p2 - p1;
    v2 = p3 - p1;
    crossz = v1(1)*v2(2) - v1(2)*v2(1);
    area = 0.5 * abs(crossz);
    % 将三角形面积均分给三个顶点（lumped mass）
    vertexArea(tri) = vertexArea(tri) + area/3;
end

% 构造对角稀疏矩阵
M = sparse(1:N, 1:N, vertexArea, N, N);
end