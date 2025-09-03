function L = cotangent_Graph_Laplacian(pts, TR) 
% COTANGENT_GRAPH_LAPLACIAN 构造余切图拉普拉斯矩阵 L_C 
%   L = cotangent_Graph_Laplacian(pts, TR)
%   输入：
%       pts - N×2 顶点坐标 [x, y]
%       TR  - delaunayTriangulation 对象
%   输出：
%       L   - N×N 稀疏余切拉普拉斯矩阵

    % 顶点数量 N 和三角形列表 T (M×3) 
    N = size(pts,1);
    % 获取所有边和对应邻接三角形索引（一次性调用提升性能）
    E = TR.edges;                         % L×2，每行是一条边的两个顶点
    adj = TR.edgeAttachments(E(:,1), E(:,2));

    % 预估非零条目数: 每个三角形3条边，每条边2个非对角 + 2个对角
    numEdges = size(E,1);
    I = zeros(numEdges*4,1);
    J = zeros(numEdges*4,1);
    V = zeros(numEdges*4,1);
    idx = 0;

    % 遍历每条边，仅处理一次
    for k = 1:numEdges
        a = E(k,1);
        b = E(k,2);
        tris = adj{k};
        if isempty(tris)
            continue;
        end
        % 第一三角形中对顶点 c
        t1 = tris(1);
        tri1 = TR.ConnectivityList(t1,:);
        c = tri1(~ismember(tri1, [a,b]));
        u = pts(a,:) - pts(c,:);
        v = pts(b,:) - pts(c,:);
        cross1 = cross([u,0],[v,0]);
        area2_1 = norm(cross1);
        if area2_1 < eps
            cot1 = 0;
        else
            cot1 = dot(u,v) / area2_1;
        end

        % 第二三角形（若存在）对顶点 d
        if numel(tris) > 1
            t2 = tris(2);
            tri2 = TR.ConnectivityList(t2,:);
            d = tri2(~ismember(tri2, [a,b]));
            % 计算 cot_beta 
            u2 = pts(a,:) - pts(d,:);
            v2 = pts(b,:) - pts(d,:);
            cross2 = cross([u2,0],[v2,0]);
            area2_2 = norm(cross2);
            if area2_2 < eps
                cot2 = 0;  % 第二三角形退化
            else
                cot2 = dot(u2,v2) / area2_2;
            end
        else
            cot2 = 0;
        end
        w = 0.5 * (cot1 + cot2);

        % 累加到稀疏矩阵条目
        idx = idx + 1; I(idx)=a; J(idx)=b; V(idx)=-w;
        idx = idx + 1; I(idx)=b; J(idx)=a; V(idx)=-w;
    end

    % 构造稀疏非对角矩阵并补对角使每行和为0
    L_off = sparse(I(1:idx), J(1:idx), V(1:idx), N, N);
    L = L_off;
    diagVals = -sum(L_off,2);
    L = L + sparse(1:N, 1:N, diagVals, N, N);
    
    % 验证
    if norm(L - L', 'fro') > 1e-8
        warning('L 非对称（超出数值容差）');
    end
end