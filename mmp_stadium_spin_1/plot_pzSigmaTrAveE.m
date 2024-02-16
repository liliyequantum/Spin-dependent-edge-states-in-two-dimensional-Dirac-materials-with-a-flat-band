clear
close all
load('data2Plot_sigmaTrAveE_a.mat')
load('data2Plot_sigmaTrAveEChaotic_a.mat')

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

% f = figure();
% f.Position = [100 100 600 400];
% plot(aa, pzAveE,'color',color_yellow,'LineWidth',2)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
% xlabel('$a$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% ylabel('$\langle P_{z} \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20);
% axis tight
% ylim([0,0.3])

f = figure();
f.Position = [100 100 600 400];
aa = 0:0.05:1.5;
plot(aa, REdgeAveE,'color',color_blue,'LineWidth',2);hold on;
scatter(aa,REdgeAveE,10,'MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', color_blue,...
              'LineWidth',1.5)
hold on
Idx = 1:1:length(aa);
plot(aa, RChaoticAveE(1,(Idx-1)*5+1),'color',color_red,'LineWidth',2);hold on;
scatter(aa,RChaoticAveE(1,(Idx-1)*5+1),10,'MarkerEdgeColor', color_red,...
              'MarkerFaceColor', color_red,...
              'LineWidth',1.5)
hold on
% plot(aaChaotic, RChaoticAveE,'color',color_red,'LineWidth',2);hold on;
% scatter(aaChaotic,RChaoticAveE,10,'MarkerEdgeColor', color_red,...
%               'MarkerFaceColor', color_red,...
%               'LineWidth',1.5)
% hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$a$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle \sigma_{tr} \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25);
axis tight
