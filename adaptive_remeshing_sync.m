function [pts_new, U_new, TR_new] = adaptive_remeshing_sync(pts, TR, U, theta_max, theta_min, K_dB)
% A simple remeshing synchronization: split edges by curvature and remove low-power verts,
% and interpolate U onto new pts using scatteredInterpolant.
    % compute approximate normals (centroid method)
    T = TR.ConnectivityList;
    N0 = size(pts,1);
    Vnorm = zeros(N0,2);
    for i = 1:N0
        triIdx = any(T==i,2);
        tris = T(triIdx,:);
        for j = 1:size(tris,1)
            tri = tris(j,:);
            centroid = mean(pts(tri,:),1);
            v1 = centroid - pts(i,:);
            n1 = [-v1(2), v1(1)];
            Vnorm(i,:) = Vnorm(i,:) + n1;
        end
        if norm(Vnorm(i,:))>0, Vnorm(i,:) = Vnorm(i,:) / norm(Vnorm(i,:)); end
    end

    newPts = pts;
    newU = U;

    E = TR.edges;
    toDelete = false(size(newPts,1),1);
    addPts = []; addU = [];

    for e = 1:size(E,1)
        a = E(e,1); b = E(e,2);
        if a>size(Vnorm,1) || b>size(Vnorm,1), continue; end
        cosTh = dot(Vnorm(a,:), Vnorm(b,:));
        cosTh = max(min(cosTh,1),-1);
        theta = acosd(cosTh);
        if theta > theta_max
            mid = mean(pts([a,b],:),1);
            midU = mean(U([a,b]));
            addPts = [addPts; mid];
            addU = [addU; midU];
        elseif theta < theta_min
            toDelete(b) = true;
        end
    end

    if ~isempty(addPts)
        newPts = [newPts; addPts];
        newU = [newU; addU];
    end

    Pd = abs(newU).^2;
    Pmean = mean(Pd);
    Kth = 10^(K_dB/10);
    mask_power = Pd >= Kth * Pmean;

    keepMask = (~toDelete) & mask_power;
    if ~any(keepMask)
        keepMask = true(size(keepMask));
    end

    pts_new = newPts(keepMask, :);
    U_new_raw = newU(keepMask);

    % interpolate (use scatteredInterpolant to transfer old U to new pts)
    F = scatteredInterpolant(pts(:,1), pts(:,2), U, 'natural', 'nearest');
    U_new = F(pts_new(:,1), pts_new(:,2));
    TR_new = delaunayTriangulation(pts_new);
end

function epsilon = compute_error_KM(pts, TR, U)
    [K, A_vert] = cotangent_K_and_M(pts, TR);
    M = spdiags(A_vert,0,size(K,1),size(K,1));
    LU = M \ (K * U);
    E = TR.edges;
    d = pts(E(:,1),:) - pts(E(:,2),:);
    len = sqrt(sum(d.^2,2));
    N = size(pts,1);
    sum_len = accumarray([E(:,1);E(:,2)], [len;len], [N,1]);
    deg = accumarray([E(:,1);E(:,2)], 1, [N,1]);
    h_i = sum_len ./ max(deg,1);
    epsilon = mean( (h_i.^2) .* abs(LU) );
end