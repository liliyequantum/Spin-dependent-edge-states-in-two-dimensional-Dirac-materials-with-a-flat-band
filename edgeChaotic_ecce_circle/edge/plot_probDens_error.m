clear
close all

% load('dataProbDensEdgeState.mat')
load('./data/dataProbDensEdgeStateAndError.mat')

figure()
pcolor(xx, yy, probDens); hold on;
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
% title(['S Matrix E=',num2str(E)],'interpreter','latex')
shading flat; 
% colorbar;
set(gca,'fontsize',25,'FontName','Times New Roman')

figure()
pcolor(xx, yy, err); hold on;
% colormap(jet(100));
plot(real(Zj1), imag(Zj1), 'w-');hold on
plot(real(Zj2), imag(Zj2), 'w-');hold on
shading flat; 
colorbar;
set(gca,'fontsize',15)