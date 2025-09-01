function [pts, TR] = generate_Fibonacci_mesh(N, R)
    phi = (1 + sqrt(5)) / 2;
    goldenAngle = 2*pi*(1 - 1/phi);
    pts = zeros(N,2);
    for k = 0:(N-1)
        r = R * sqrt((k + 0.5) / N);
        theta = k * goldenAngle;
        pts(k+1,1) = r * cos(theta);
        pts(k+1,2) = r * sin(theta);
    end
    TR = delaunayTriangulation(pts(:,1), pts(:,2));
end