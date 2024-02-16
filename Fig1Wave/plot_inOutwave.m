clear
close all

mu = 24;
xi = 0;
E = 15.075;
Lrange = 35;

N = 600;
xyra = 1.5;

% wave function
load(['./data/data2plotInOutWave_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'])

figure()
pcolor(xx, yy, probIn); hold on;
title(['S Matrix E=',num2str(E),' probIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

f = figure();
f.Position = [10,10,300,250];
pcolor(xx, yy, probOut); hold on;
%title(['S Matrix E=',num2str(E),' probScatter'],'interpreter','latex')
shading flat; 
colorbar;
clim([0 1])
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, prob); hold on;
title(['S Matrix E=',num2str(E),' prob'],'interpreter','latex')
shading flat; 
colorbar;
% clim([0 10])
set(gca,'fontsize',25,'FontName','times new roman')