function U1 = propagate_HalfStep_mass(U0, Psi, lambda_vec, M, dz, k)
    U_hat = Psi' * (M * U0);                     % M-weighted projection
    % Phase factor per paper: exp(j * lambda * dz/(2*k))
    phase = exp(1i * lambda_vec * (dz/(2*k)));
    U_hat = phase .* U_hat;
    U1 = Psi * U_hat;
end