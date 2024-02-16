clear
close all

dbstop if warning
% dbclear if warning

xi = 0.12;
spinIndex = 1; % up = 1, down = -1

if spinIndex==-1
    v1 = -10; 
    v2 = 40;
    L = -40:1:40;
    EE = 14.5:1e-4/2:15.3;
else
    mu = 22;
    v1 = -10 - 2*mu;
    v2 = 40 - 2*mu;
    L = -35:1:35;
    EE = 14.5:1e-3/2:15.3;
end

R_1 = 1;
R_2 = 0.6;

A_0 = diag(0.5/sqrt(2)*1j.^(L-1));
sigma_bar_total = zeros(1,length(EE));
sigma_total = zeros(1,length(EE));
for NN=1:1:length(EE)
    E = EE(NN);

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

    sigma_bar_total(NN) = 0.5*sum(sum(abs(S - eye(size(S))).^2))./k_0;
    TM = S - eye(size(S));
    sigma_total(NN) = transpose(diag(A_0))*(TM*TM')*conj(diag(A_0))*4/k_0;
    disp(NN/length(EE))
end

save Fig1_chaotic_sigmaBarTotal_mu_22_supp


disp(['imag of sigma_bar_t = ',num2str(max(abs(imag(sigma_bar_total))))])
disp(['imag of sigma_t = ',num2str(max(abs(imag(sigma_total))))])

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

f = figure(1);
f.Position = [100 100 600 400];

if spinIndex==1 % spin up
    plot(EE,real(sigma_bar_total),'color',color_red,'LineWidth',3)
    hold on
else % spin down
    plot(EE,real(sigma_bar_total),'color',color_blue,'LineWidth',3)
    hold on 
end
% ylim([1.9,2.05])
% title(['$\xi = $',num2str(xi)],'Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
xlabel('$E$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
ylabel('$\bar{\sigma}_t$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
axis tight

f = figure(2);
f.Position = [100 100 600 400];
if spinIndex==1 % spin up
    plot(EE,real(sigma_total),'color',color_red,'LineWidth',2)
    hold on
else % spin down
    plot(EE,real(sigma_total),'color',color_blue,'LineWidth',2)
    hold on
end
title(['$\xi = $',num2str(xi)],'Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
xlabel('$E$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
ylabel('$\sigma_t$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
axis tight

figure()
pcolor(sqrt(sqrt(abs(S)))); hold on;
colormap(jet(20));
shading flat; 
colorbar;
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)

numL = (length(L)-1)/2;
SArray = diag(S);
Sx = SArray(1:1:numL);
Sy = SArray(end:-1:(numL+2));
figure()
plot([-1,1],[-1,1],'k-','linewidth',1);hold on;
scatter(real(Sx),real(Sy),100,'MarkerEdgeColor', color_red,...
              'MarkerFaceColor', ["none"],...
              'LineWidth',1.5);hold on;
scatter(imag(Sx),imag(Sy),100,'MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', ["none"],...
              'LineWidth',1.5);hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$S_{-l,-l}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$S_{l,l}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
