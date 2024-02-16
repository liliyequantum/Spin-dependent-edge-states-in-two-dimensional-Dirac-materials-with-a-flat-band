clear
close all
% maxNumCompThreads(15)

xi = 0.03;
E = 15.1;

mu = 24;
v1 = -10 - 2*mu;
v2 = 40 - 2*mu;
L = -35:1:35;

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
    num2str(xi),'_E_',num2str(E),'.mat'])

% %% boundary check
% 
% % boundary1
% 
% Nj1=400;
% theta1=linspace(0,2*pi,Nj1+1);
% theta1=theta1(1:Nj1);
% Zj1=R_1*exp(sqrt(-1)*theta1);
% 
% % boundary2
% 
% Nj2=400;
% theta2op=linspace(0,2*pi,Nj2+1);
% theta2op=theta2op(1:Nj2);
% Zj2op=R_2*exp(sqrt(-1)*theta2op);
% Zj2 = xi + Zj2op;
% 
% r1 = abs(Zj1);
% theta1 = angle(Zj1);
% % r1p = abs(Zj1-xi);
% % theta1p = angle(Zj1-xi);
% r2p = abs(Zj2op);
% theta2p = angle(Zj2op);
% 
% % psi_I at boundary 1
% % disp('psi_I, o origin, at boundary 1')
% psi_I_1=zeros(1,Nj1);
% psi_I_2=zeros(1,Nj1);
% psi_I_3=zeros(1,Nj1);
% for mm = 1:1:length(L)
%         for ll = 1:1:length(L)
%             if ll==mm
%                 psi_I_1 = psi_I_1 + A_0(mm,mm)*besselh(L(mm)-1, 2, k_0*r1).*exp(1j*(L(mm)-1)*theta1);
%                 psi_I_2 = psi_I_2 + A_0(mm,mm)*1j*alpha_0*sqrt(2)*besselh(L(mm), 2, k_0*r1).*exp(1j*L(mm)*theta1);
%                 psi_I_3 = psi_I_3 - A_0(mm,mm)*besselh(L(mm)+1, 2, k_0*r1).*exp(1j*(L(mm)+1)*theta1);
%             end
%             psi_I_1 = psi_I_1 + A_0(mm,mm)*S(mm,ll)*besselh(L(ll)-1, 1, k_0*r1).*exp(1j*(L(ll)-1)*theta1);
%             psi_I_2 = psi_I_2 + A_0(mm,mm)*S(mm,ll)*1j*alpha_0*sqrt(2)*besselh(L(ll), 1, k_0*r1).*exp(1j*L(ll)*theta1);
%             psi_I_3 = psi_I_3 - A_0(mm,mm)*S(mm,ll)*besselh(L(ll)+1, 1, k_0*r1).*exp(1j*(L(ll)+1)*theta1);
%         end
%         disp({'psi_I at boundary 1: ', num2str(mm/length(L))})
% end
% 
% % psi_II at boundary 1
% % disp('psi_II, o origin, at boundary 1')
% psi1_II_1=zeros(1,Nj1);
% psi1_II_2=zeros(1,Nj1);
% psi1_II_3=zeros(1,Nj1);
% for mm = 1:1:length(L)
%     for ll = 1:1:length(L)
%         psi1_II_1 = psi1_II_1 + A(mm,ll)*besselh(L(ll)-1, 2, k_1*r1).*exp(1j*(L(ll)-1)*theta1);
%         psi1_II_2 = psi1_II_2 + A(mm,ll)*1j*alpha_1*sqrt(2)*besselh(L(ll), 2, k_1*r1).*exp(1j*L(ll)*theta1);
%         psi1_II_3 = psi1_II_3 - A(mm,ll)*besselh(L(ll)+1, 2, k_1*r1).*exp(1j*(L(ll)+1)*theta1);
%         for llp = 1:1:length(L)
%             psi1_II_1 = psi1_II_1 + A(mm,ll)*S_od(ll,llp)*besselh(L(llp)-1, 1, k_1*r1).*exp(1j*(L(llp)-1)*theta1);
%             psi1_II_2 = psi1_II_2 + A(mm,ll)*S_od(ll,llp)*1j*alpha_1*sqrt(2)*besselh(L(llp), 1, k_1*r1).*exp(1j*L(llp)*theta1);
%             psi1_II_3 = psi1_II_3 - A(mm,ll)*S_od(ll,llp)*besselh(L(llp)+1, 1, k_1*r1).*exp(1j*(L(llp)+1)*theta1);
%         end
%     end
%     disp({'psi_II at boundary 1: ', num2str(mm/length(L))})
% end
% 
% % psi_II at boundary 2
% disp('psi_II, oprime origin, at boundary 1')
% psi2_II_1=zeros(1,Nj2);
% psi2_II_2=zeros(1,Nj2);
% psi2_II_3=zeros(1,Nj2);
% for mm = 1:1:length(L)
%     for ll = 1:1:length(L)
% 
%         psi2_II_1 = psi2_II_1 + A_tilde(mm,ll)*(besselh(L(ll)-1, 2, k_1*r2p).*exp(1j*(L(ll)-1)*theta2p)+...
%             S_cd(ll,ll)*besselh(L(ll)-1, 1, k_1*r2p).*exp(1j*(L(ll)-1)*theta2p));
%         psi2_II_2 = psi2_II_2 + A_tilde(mm,ll)*1j*alpha_1*sqrt(2)*(besselh(L(ll), 2, k_1*r2p).*exp(1j*L(ll)*theta2p)+...
%             S_cd(ll,ll)*besselh(L(ll), 1, k_1*r2p).*exp(1j*L(ll)*theta2p));
%         psi2_II_3 = psi2_II_3 - A_tilde(mm,ll)*(besselh(L(ll)+1, 2, k_1*r2p).*exp(1j*(L(ll)+1)*theta2p)+...
%             S_cd(ll,ll)*besselh(L(ll)+1, 1, k_1*r2p).*exp(1j*(L(ll)+1)*theta2p));
% 
%     end
%     disp({'psi_II at boudnary 2: ', num2str(mm/length(L))})
% end
% 
% % psi_III
% disp('psi_III, oprime origin')
% psi_III_1=zeros(1,Nj2);
% psi_III_2=zeros(1,Nj2);
% psi_III_3=zeros(1,Nj2);
% for mm = 1:1:length(L)
%     for ll = 1:1:length(L)
%         psi_III_1 = psi_III_1 + B_tilde(mm,ll)*besselj(L(ll)-1, k_2*r2p).*exp(1j*(L(ll)-1)*theta2p);
%         psi_III_2 = psi_III_2 + B_tilde(mm,ll)*1j*alpha_2*sqrt(2)*besselj(L(ll), k_2*r2p).*exp(1j*L(ll)*theta2p);
%         psi_III_3 = psi_III_3 - B_tilde(mm,ll)*besselj(L(ll)+1, k_2*r2p).*exp(1j*(L(ll)+1)*theta2p);
%     end
%     disp({'psi_III at boundary 3: ', num2str(mm/length(L))})
% end
% 
% psi_I_1=psi_I_1/sqrt(2);
% psi_I_2=psi_I_2/sqrt(2);
% psi_I_3=psi_I_3/sqrt(2);
% psi1_II_1=psi1_II_1/sqrt(2);
% psi1_II_2=psi1_II_2/sqrt(2);
% psi1_II_3=psi1_II_3/sqrt(2);
% psi2_II_1=psi2_II_1/sqrt(2);
% psi2_II_2=psi2_II_2/sqrt(2);
% psi2_II_3=psi2_II_3/sqrt(2);
% psi_III_1=psi_III_1/sqrt(2);
% psi_III_2=psi_III_2/sqrt(2);
% psi_III_3=psi_III_3/sqrt(2);
% 
% left1 = psi_I_2;
% right1 = psi1_II_2;
% left2 = psi_I_1.*exp(1i*theta1) + psi_I_3.*exp(-1i*theta1);
% right2 = psi1_II_1.*exp(1i*theta1) + psi1_II_3.*exp(-1i*theta1);
% error_1 = (abs(left1-right1)+abs(left2-right2))./(abs(left1)+abs(left2)+abs(right1)+abs(right2));
% 
% left1 = psi_III_2;
% right1 = psi2_II_2;
% left2 = psi_III_1.*exp(1i*theta2op) + psi_III_3.*exp(-1i*theta2op);
% right2 = psi2_II_1.*exp(1i*theta2op) + psi2_II_3.*exp(-1i*theta2op);
% error_2 = (abs(left1-right1)+abs(left2-right2))./(abs(left1)+abs(left2)+abs(right1)+abs(right2));
% 
% figure()
% semilogy(error_1);hold on;
% semilogy(error_2);hold on;
% set(gca,'fontsize',25)
% xlabel('boundary')
% legend('boundary 1', 'boundary 2')

