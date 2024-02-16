clear
close all

mu = 22;
xi = 0.12;
E = 14.8;
Lrange = 35;

N = 450;
xyra = 1.1;

% wave function
load(['./data/data2plotDensityAndError_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'])

f = figure();
f.Position=[10,10,300,250];
pcolor(xx, yy, probDens); hold on;
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
%title(['S Matrix E=',num2str(E)],'interpreter','latex')
shading flat; 
%colorbar;
%set(gca,'fontsize',35,'FontName','times new roman')
axis off

figure()
pcolor(xx, yy, err); hold on;
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
shading flat; 
colorbar;
set(gca,'fontsize',15)