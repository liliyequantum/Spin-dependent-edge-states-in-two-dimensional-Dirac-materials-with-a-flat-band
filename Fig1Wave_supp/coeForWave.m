function coeForWave(mu,xi,E,Lrange)

    % xi = 0.03;
    % E = 15.1;
    % 
    % mu = 24;
    v1 = -10 - 2*mu;
    v2 = 40 - 2*mu;
    L = -Lrange:1:Lrange;
    
    R_1 = 1;
    R_2 = 0.6;
    
    k_0 = abs(E);
    k_1 = abs(E - v1);
    k_2 = abs(E - v2);
    alpha_0 = sign(E);
    alpha_1 = sign(E - v1);
    alpha_2 = sign(E - v2);
    
    A_0 = diag(0.5/sqrt(2)*1j.^(L-1));
    
    X_1 = diag(besselh(L, 1, k_0*R_1));
    X_2 = diag(besselh(L, 2, k_0*R_1));
    x_1 = diag(besselh(L, 1, k_1*R_1));
    x_2 = diag(besselh(L, 2, k_1*R_1));
    
    Y_1 = diag(besselh(L+1, 1, k_0*R_1));
    Y_2 = diag(besselh(L+1, 2, k_0*R_1));
    y_1 = diag(besselh(L+1, 1, k_1*R_1));
    y_2 = diag(besselh(L+1, 2, k_1*R_1));
    
    Z_1 = diag(besselh(L-1, 1, k_0*R_1));
    Z_2 = diag(besselh(L-1, 2, k_0*R_1));
    z_1 = diag(besselh(L-1, 1, k_1*R_1));
    z_2 = diag(besselh(L-1, 2, k_1*R_1));
    
    temp_1 = besselj(L, k_2*R_2).*(besselh(L-1, 2, k_1*R_2) -besselh(L+1, 2, k_1*R_2)) ...
             - alpha_1*alpha_2*besselh(L, 2, k_1*R_2).*(besselj(L-1, k_2*R_2) - besselj(L+1, k_2*R_2));
    temp_2 =   besselj(L, k_2*R_2).*(besselh(L-1, 1, k_1*R_2) -besselh(L+1, 1, k_1*R_2)) ...
             - alpha_1*alpha_2*besselh(L, 1, k_1*R_2).*(besselj(L-1, k_2*R_2) - besselj(L+1, k_2*R_2));
    S_cd = diag(-temp_1./temp_2);
    [column,row] = meshgrid(L,L);
    U = besselj(column - row, k_1*xi);
    U_inv = besselj(row - column, k_1*xi);
    S_od = U_inv*S_cd*U;
    
    F = x_2 + S_od*x_1;
    G = y_2 + S_od*y_1;
    H = z_2 + S_od*z_1;
    T = F\(H-G);
    
    S = -(Z_2 - Y_2-alpha_0*alpha_1*X_2*T)/(Z_1 - Y_1 - alpha_0*alpha_1*X_1*T);
    A = alpha_0*alpha_1*A_0*(X_2+S*X_1)/F;
    
    J = besselj(row - column, k_1*xi);
    A_tilde = A*J;
    B_tilde = alpha_1*alpha_2*A_tilde*diag((besselh(L, 2, k_1*R_2)...
        +transpose(diag(S_cd)).*besselh(L, 1, k_1*R_2))./besselj(L, k_2*R_2));
    
    save(['./data/coeDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'.mat'])
end

