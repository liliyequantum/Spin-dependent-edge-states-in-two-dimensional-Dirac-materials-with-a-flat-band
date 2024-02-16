clear
close all
load('./pictureData/dataToPlot_sigmaTr_G_Pz_mu17p5to27p5MoreFine_xi.mat')

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

if exist('picture','dir')==0
    mkdir('picture')
end

% plot
[XiM,muM] = meshgrid(Xi, mmu); 
figure()
pcolor(XiM, muM, pzAveE); hold on;
%title('$\langle P_{z} \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,'./picture/pz.pdf')

f = figure();
% f.Position = [10,10,600,300];
plot(Xi,pzAveE(39,:),'color',color_yellow,'LineWidth',3)
hold on
scatter(Xi,pzAveE(39,:),10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',2.5)
hold on
% title(['$\mu=$',num2str(mmu(39))],'Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle P_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
xticks([0,0.1,0.2,0.3])
yticks([0.5,0.7,0.97])
% ylim([0,1])
grid on

f = figure();
plot(mmu,pzAveE(:,1),'color',color_yellow,'LineWidth',3)
hold on
scatter(mmu,pzAveE(:,1),10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',2.5)
hold on
% title('$\xi=0$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',35)
set(gca, 'LineWidth',3,'Fontname', 'Times New Roman','FontSize',45)
% xlabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',35)
% ylabel('$\langle P_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',35)
axis tight
ylim([0,1])
xticks([20,24])
% saveas(gcf,['./picture/mu_',num2str(mu),'/pz.pdf'])
grid on

maxPzAveE = max(transpose(pzAveE));
figure()
plot(mmu,maxPzAveE,'color',color_purple,'LineWidth',2)
hold on
scatter(mmu,maxPzAveE,10,'MarkerEdgeColor', color_purple,...
              'MarkerFaceColor', color_purple,...
              'LineWidth',1.5)
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\max\langle P_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
saveas(gcf,'./picture/maxPz.pdf')
% ylim([0.6,1])
axis tight
yticks([0.7,0.8,0.87,0.9,0.97])

figure()
pcolor(XiM, muM, pzGAveE); hold on;
title('$\langle P_{z} \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25)
title('based on G')
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,'./picture/pzG.pdf')

figure()
pcolor(XiM, muM, GChaoticAveE); hold on;
title('$\langle G\rangle$ Chaotic','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,'./picture/GChaotic.pdf')

figure()
pcolor(XiM, muM, sigmaTrChaoticAveE); hold on;
title('$\langle \sigma_{tr}\rangle$ Chaotic','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,'./picture/sigmaTrChaotic.pdf')

% plot(Xi,pzAveE(27,:))