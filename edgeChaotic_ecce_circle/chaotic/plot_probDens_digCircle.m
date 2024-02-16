clear
close all

mu = 22;
xi = 0.12;
E = 14.8;
 
N = 500;
xyra = 3;

Rj=1.5;
Nj=300;
theta=linspace(0,2*pi,Nj+1);
theta=theta(1:Nj);
Zj=Rj*exp(sqrt(-1)*theta);

x_choose=linspace(-xyra,xyra,N);
y_choose=linspace(-xyra,xyra,N);
[xx,yy]=meshgrid(x_choose,y_choose);
zz=xx+sqrt(-1)*yy;

Out=ones(N,N);
for i=1:N
    for j=1:N
        if inpolygon(xx(i,j),yy(i,j),real(Zj),imag(Zj))   
            Out(i,j)=0;
        end
    end
    disp(i/N)
end

% wave function
load(['./data/data2plotDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'])

figure()
pcolor(xx, yy, probDens.*Out); hold on;
% colormap(jet(20));
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
title(['S Matrix E=',num2str(E)],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')
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
% figure()
% pcolor(xx, yy, err); hold on;
% % colormap(jet(20));
% plot(real(Zj1), imag(Zj1), 'w-');hold on
% plot(real(Zj2), imag(Zj2), 'w-');hold on
% shading flat; 
% colorbar;
% set(gca,'fontsize',15)