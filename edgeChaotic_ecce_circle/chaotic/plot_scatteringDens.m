clear
close all

mu = 24;
xi = 0;
E = 14.8;

N=10;
xleft = 1;%1.1;
xright = 6;

% wave function
load(['./data/data2plotDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_N_',num2str(N),'_xleft_',...
        num2str(xleft),'_xright_',num2str(xright),'.mat'])

figure()
pcolor(xx, yy, probDens); hold on;
title(['S Matrix E=',num2str(E)],'interpreter','latex')
shading flat; 
% colorbar;
set(gca,'fontsize',25,'FontName','times new roman')
