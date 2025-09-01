function epsilon = compute_error_KM(pts, TR, U)
% *** MODIFIED: compute error using K and M (discrete Laplacian M\K*U) ***
    [K, A_vert] = cotangent_K_and_M(pts, TR);
    M = spdiags(A_vert,0,size(K,1),size(K,1));
    LU = M \ (K * U);
    % estimate average edge length per vertex
    E = TR.edges;
    d = pts(E(:,1),:) - pts(E(:,2),:);
    len = sqrt(sum(d.^2,2));
    N = size(pts,1);
    sum_len = accumarray([E(:,1);E(:,2)], [len;len], [N,1]);
    deg = accumarray([E(:,1);E(:,2)], 1, [N,1]);
    h_i = sum_len ./ max(deg,1);
    epsilon = mean( (h_i.^2) .* abs(LU) );
end