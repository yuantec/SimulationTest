function [K, A_vert] = cotangent_K_and_M(pts, TR)
% returns stiffness K and lumped mass A_vert
    T = TR.ConnectivityList;
    N = size(pts,1);
    E = TR.edges;
    adj = TR.edgeAttachments(E(:,1), E(:,2));

    I=[]; J=[]; V=[];
    for kk = 1:size(E,1)
        a = E(kk,1); b = E(kk,2);
        tris = adj{kk};
        if isempty(tris), continue; end

        % first triangle
        t1 = tris(1); tri1 = T(t1,:);
        c = tri1(~ismember(tri1,[a,b]));
        u = pts(a,:) - pts(c,:); v = pts(b,:) - pts(c,:);
        area2_1 = norm(cross([u,0],[v,0]));
        cot1 = 0; if area2_1>eps, cot1 = dot(u,v)/area2_1; end

        % second triangle if exists
        if numel(tris) > 1
            t2 = tris(2); tri2 = T(t2,:);
            d = tri2(~ismember(tri2,[a,b]));
            u2 = pts(a,:) - pts(d,:); v2 = pts(b,:) - pts(d,:);
            area2_2 = norm(cross([u2,0],[v2,0]));
            cot2 = 0; if area2_2>eps, cot2 = dot(u2,v2)/area2_2; end
        else
            cot2 = 0;
        end

        w = 0.5 * (cot1 + cot2);   % keep 1/2 factor

        I(end+1)=a; J(end+1)=b; V(end+1)=-w;
        I(end+1)=b; J(end+1)=a; V(end+1)=-w;
        I(end+1)=a; J(end+1)=a; V(end+1)= w;
        I(end+1)=b; J(end+1)=b; V(end+1)= w;
    end

    K = sparse(I,J,V,N,N);
    K = (K + K')/2;

    % lumped mass: each triangle area / 3
    areas = zeros(size(T,1),1);
    for t = 1:size(T,1)
        ids = T(t,:);
        p1 = pts(ids(1),:); p2 = pts(ids(2),:); p3 = pts(ids(3),:);
        areas(t) = 0.5 * abs( (p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)) );
    end
    A_vert = accumarray(T(:), repmat(areas/3,3,1), [N,1]);
    A_vert(A_vert<=0) = eps;
end