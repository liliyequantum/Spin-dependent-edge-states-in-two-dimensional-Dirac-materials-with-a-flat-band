clear
close all
mu = 24;
load(['./pictureData/dataToPlot_sigmaTrOverTotal_Pz_mu_',num2str(mu),'.mat'])

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
plot(Xi,sigmaTrTotalEdgeAveE,'color',color_blue,'LineWidth',2)
hold on
scatter(Xi,sigmaTrTotalEdgeAveE,10,'MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', color_blue,...
              'LineWidth',1.5)
hold on
plot(Xi,sigmaTrTotalChaoticAveE,'color',color_red,'LineWidth',2)
hold on
scatter(Xi,sigmaTrTotalChaoticAveE,10,'MarkerEdgeColor', color_red,...
              'MarkerFaceColor', color_red,...
              'LineWidth',1.5)
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle\sigma_{tr}/\sigma_{t}\rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight


f = figure();
% f.Position = [100 100 500 500];
plot(Xi,pzOverTotalAveE,'color',color_yellow,'LineWidth',2)
hold on
scatter(Xi,pzOverTotalAveE,10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',1.5)
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle P_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
ylim([0.4,1])
