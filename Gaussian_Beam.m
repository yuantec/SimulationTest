function U0 = Gaussian_Beam(pts, w0, z, k)
  % 输入：
  %   pts - N×2 网格顶点坐标 [x,y]
  %   w0  - 复束腰半径
  %   z   - 初始化平面位置 (m)
  %   k   - 波数 = 2*pi/lambda
  %
  % 输出：
  %   U0  - N×1 复振幅初值
  
  % 1. 计算瑞利距离
  ZR = k * w0^2 / 2;
  
  % 2. 束宽 w(z)
  wz = w0 * sqrt(1 + (z/ZR)^2);
  
  % 3. 曲率半径 R(z)
  if z == 0
    Rz = Inf;    % 平面波前
  else
    Rz = z * (1 + (ZR/z)^2);
  end
  
  % 4. Gouy 相位
  PhiG = atan(z/ZR);
  
  % 5. 计算每个顶点的复振幅
  x = pts(:,1);
  y = pts(:,2);
  r2 = x.^2 + y.^2;
  
  % 振幅衰减因子
  amp = (w0 / wz) * exp( -r2 / (wz^2) );
  
  % 二次相位项；平面波前时该项为 0
  if isfinite(Rz)
    quadPhase = exp(-1i * k * r2 / (2*Rz));
  else
    quadPhase = 1;
  end
  
  % Gouy 相位
  gouyPhase = exp(-1i * PhiG);
  
  % 合成初始场
  U0 = amp .* quadPhase .* gouyPhase;
end
