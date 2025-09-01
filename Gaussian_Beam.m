function U0 = Gaussian_Beam(pts, w0, z, k)
    ZR = k * w0^2 / 2;
    wz = w0 * sqrt(1 + (z/ZR)^2);
    if z==0, Rz = Inf; else Rz = z * (1 + (ZR/z)^2); end
    PhiG = atan(z/ZR);
    x = pts(:,1); y = pts(:,2);
    r2 = x.^2 + y.^2;
    amp = (w0 / wz) * exp( -r2 / (wz^2) );
    if isfinite(Rz)
        quadPhase = exp(-1i * k * r2 / (2*Rz));
    else
        quadPhase = 1;
    end
    gouyPhase = exp(1i * PhiG);
    U0 = amp .* quadPhase .* gouyPhase;
end