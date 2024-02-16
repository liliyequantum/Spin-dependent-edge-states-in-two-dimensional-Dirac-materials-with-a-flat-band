clear
close all
mu = 20;
load(['./pictureData/dataToPlot_sigmaTr_G_Pz_mu_',num2str(mu),'.mat'])

% plot

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

if exist(['./picture/mu_',num2str(mu)],'dir')==0
    mkdir(['./picture/mu_',num2str(mu)]);
end

figure()
plot(Xi,sigmaTrEdgeAveE,'color',color_blue,'LineWidth',2)
hold on
scatter(Xi,sigmaTrEdgeAveE,10,'MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', color_blue,...
              'LineWidth',1.5)
hold on
plot(Xi,sigmaTrChaoticAveE,'color',color_red,'LineWidth',2)
hold on
scatter(Xi,sigmaTrChaoticAveE,10,'MarkerEdgeColor', color_red,...
              'MarkerFaceColor', color_red,...
              'LineWidth',1.5)
hold on
yticks([0.4,0.8,1.2])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle\sigma_{tr}\rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,['./picture/mu_',num2str(mu),'/sigmaTr.pdf'])

figure()
plot(Xi,GEdgeAveE,'color',color_blue,'LineWidth',2)
hold on
scatter(Xi,GEdgeAveE,10,'MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', color_blue,...
              'LineWidth',1.5)
hold on
plot(Xi,GChaoticAveE,'color',color_red,'LineWidth',2)
hold on
scatter(Xi,GChaoticAveE,10,'MarkerEdgeColor', color_red,...
              'MarkerFaceColor', color_red,...
              'LineWidth',1.5)
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle G \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,['./picture/mu_',num2str(mu),'/G.pdf'])

f = figure();
% f.Position = [100 100 500 500];
plot(Xi,pzAveE,'color',color_yellow,'LineWidth',2)
hold on
scatter(Xi,pzAveE,10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',1.5)
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle P_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
ylim([0.4,1])
saveas(gcf,['./picture/mu_',num2str(mu),'/pz.pdf'])

figure()
plot(Xi,pzGAveE,'color',color_yellow,'LineWidth',2)
hold on
scatter(Xi,pzGAveE,10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',1.5)
hold on
title('based on G')
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle P_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
saveas(gcf,['./picture/mu_',num2str(mu),'/pzG.pdf'])