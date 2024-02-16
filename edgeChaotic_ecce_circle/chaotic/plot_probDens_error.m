clear
close all

mu = 24;
xi = 0;
E = 15.083;
 
N = 1200;
xyra = 3;

% wave function
load(['./data/data2plotDensityAndError_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'])
% generateIn123(xi,N,xyra)
% load(['In123_xi',num2str(xi),'_N',num2str(N),'_xyra',num2str(xyra),'.mat'])

% nright = 100;
% figure()
% pcolor(xx(:,1:nright), yy(:,1:nright), probDens(:,1:nright)); hold on;
% % colormap(jet(20));
% plot(real(Zj1), imag(Zj1), 'w-');hold on
% plot(real(Zj2), imag(Zj2), 'w-');hold on
% title(['S Matrix E=',num2str(E)],'interpreter','latex')
% shading flat; 
% colorbar;
% set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, probDens); hold on;
% colormap(jet(20));
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
title(['S Matrix E=',num2str(E)],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

% nleft = 400;
% figure()
% pcolor(xx(:,nleft:end), yy(:,nleft:end), probDens(:,nleft:end)); hold on;
% % colormap(jet(20));
% plot(real(Zj1), imag(Zj1), 'w-');hold on
% plot(real(Zj2), imag(Zj2), 'w-');hold on
% title(['S Matrix E=',num2str(E)],'interpreter','latex')
% shading flat; 
% colorbar;
% set(gca,'fontsize',25,'FontName','times new roman')




% figure()
% pcolor(xx, yy, probDens.*In1); hold on;
% % colormap(jet(20));
% plot(real(Zj1), imag(Zj1), 'w-');hold on
% plot(real(Zj2), imag(Zj2), 'w-');hold on
% title(['S Matrix E=',num2str(E)],'interpreter','latex')
% shading flat; 
% % colorbar;
% set(gca,'fontsize',25,'FontName','times new roman')
% 
% figure()
% pcolor(xx, yy, probDens.*In2); hold on;
% % colormap(jet(20));
% plot(real(Zj1), imag(Zj1), 'w-');hold on
% plot(real(Zj2), imag(Zj2), 'w-');hold on
% title(['S Matrix E=',num2str(E)],'interpreter','latex')
% shading flat; 
% % colorbar;
% set(gca,'fontsize',25,'FontName','times new roman')
% 
% figure()
% pcolor(xx, yy, probDens.*In3); hold on;
% % colormap(jet(20));
% plot(real(Zj1), imag(Zj1), 'w-');hold on
% plot(real(Zj2), imag(Zj2), 'w-');hold on
% title(['S Matrix E=',num2str(E)],'interpreter','latex')
% shading flat; 
% % colorbar;
% set(gca,'fontsize',25,'FontName','times new roman')
% 
% % error
% load(['./data/data2plotDensityError_mu_',num2str(mu),'_xi_',...
%         num2str(xi),'_E_',num2str(E),'_N_',num2str(N),'_xyra_',...
%         num2str(xyra),'.mat'])
figure()
pcolor(xx, yy, err); hold on;
% colormap(jet(20));
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
shading flat; 
colorbar;
set(gca,'fontsize',15)