clear
close all

mu = 17:1:28;
v1 = -10 - 2*mu;
v2 = 40 - 2*mu;

E = 14.9;
n1 = (E - v1)/E;
n2 = (E - v2)/E;
n1d = (E-(-10)*ones(size(mu)))/E;
n2d = (E-(40)*ones(size(mu)))/E;

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

f = figure();
% f.Position = [10,10,400,600];
plot(mu,v1,'color',color_red,'LineWidth',3);hold on;
plot(mu,v2,'color',color_red,'LineWidth',3);hold on;
plot(mu,-10*ones(size(mu)),'color',color_blue,'LineWidth',3);hold on;
plot(mu,40*ones(size(mu)),'color',color_blue,'LineWidth',3);hold on;
scatter(mu,v1,100,'>','MarkerEdgeColor', color_red,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
scatter(mu,v2,100,'o','MarkerEdgeColor', color_red,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
scatter(mu,-10*ones(size(mu)),100,'>','MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
scatter(mu,40*ones(size(mu)),100,'o','MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
% legend('$V_{1}$','$V_{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
xlabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$V$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
grid on

f=figure();
% f.Position=[10,10,400,600];
plot(mu,n1d,'color',color_blue,'LineWidth',3);hold on;
plot(mu,n2d,'color',color_blue,'LineWidth',3);hold on;
plot(mu,n1,'color',color_red,'LineWidth',3);hold on;
plot(mu,n2,'color',color_red,'LineWidth',3);hold on;
scatter(mu,n1d,100,'>','MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
scatter(mu,n2d,100,'o','MarkerEdgeColor', color_blue,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
scatter(mu,n1,100,'>','MarkerEdgeColor', color_red,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
scatter(mu,n2,100,'o','MarkerEdgeColor', color_red,...
              'MarkerFaceColor', 'white',...
              'LineWidth',3);hold on;
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
% legend('$n_{1}$','$n_{2}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
xlabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$n$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
grid on
ylim([-2,5.5])
yticks([-2,0,1,2,4,5])
