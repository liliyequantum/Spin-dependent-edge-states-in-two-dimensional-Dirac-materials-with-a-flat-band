clear
close all

dbstop if warning
% dbclear if warning

% xi = 0.165;
Xi = 0;
spinIndex = 1; % up = 1, down = -1

if spinIndex==-1
    v1 = -10; 
    v2 = 40;
    L = -40:1:40;
    EE = 14.5:1e-4/2:15.3;
else
    mu = 20;
    v1 = -10 - 2*mu;
    v2 = 40 - 2*mu;
    L = -35:1:35;
    %EE = 14.5:1e-3/2:15.3;
    E = 15;
end

R_1 = 1;
R_2 = 0.6;

A_0 = diag(0.5/sqrt(2)*1j.^(L-1));

k_0 = abs(E);
k_1 = abs(E - v1);
k_2 = abs(E - v2);
alpha_0 = sign(E);
alpha_1 = sign(E - v1);
alpha_2 = sign(E - v2);

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

Sring = zeros(length(L),length(L));
for lIdx = 1:1:length(L)
    xl = besselh(L(lIdx), 2, k_1*R_1) + besselh(L(lIdx), 1, k_1*R_1)*S_cd(lIdx,lIdx);
    yl = (besselh(L(lIdx)-1, 2, k_1*R_1)-besselh(L(lIdx)+1, 2, k_1*R_1)) + ...
        (besselh(L(lIdx)-1, 1, k_1*R_1)-besselh(L(lIdx)+1, 1, k_1*R_1))*S_cd(lIdx,lIdx);
    temp_1 = alpha_0*yl*besselh(L(lIdx), 2, k_0*R_1) - alpha_1*xl*(besselh(L(lIdx)-1, 2, k_0*R_1) - besselh(L(lIdx)+1, 2, k_0*R_1));
    temp_2 = alpha_0*yl*besselh(L(lIdx), 1, k_0*R_1) - alpha_1*xl*(besselh(L(lIdx)-1, 1, k_0*R_1) - besselh(L(lIdx)+1, 1, k_0*R_1));
    Sring(lIdx,lIdx) = -temp_1/temp_2;
end

Smatrix = zeros(length(L),length(L),length(Xi));
Sdiff = zeros(length(Xi));
for xiIdx=1:1:length(Xi)
    xi = Xi(xiIdx);

    U = besselj(column - row, k_1*xi);
    U_inv = besselj(row - column, k_1*xi);
    S_od = U_inv*S_cd*U;
    
    F = x_2 + S_od*x_1;
    G = y_2 + S_od*y_1;
    H = z_2 + S_od*z_1;
    T = F\(H-G);
    
    Smatrix(:,:,xiIdx) = -(Z_2 - Y_2-alpha_0*alpha_1*X_2*T)/(Z_1 - Y_1 - alpha_0*alpha_1*X_1*T);
    abs(Smatrix(:,:,xiIdx) - Sring)
    disp(xiIdx/length(Xi))
end

max(max(abs(Smatrix(:,:,1) - Sring)))
