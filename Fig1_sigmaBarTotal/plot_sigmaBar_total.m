clear
close all

% load('Fig1_edge_sigmaBarTotal.mat')
load('Fig1_chaotic_sigmaBarTotal.mat')

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

f = figure(1);
f.Position = [100 100 800 400];

if spinIndex==1 % spin up
    plot(EE,real(sigma_bar_total),'color',color_red,'LineWidth',3)
    hold on
else % spin down
    plot(EE,real(sigma_bar_total),'color',color_blue,'LineWidth',3)
    hold on 
end
set(gca,'box','off')
% ylim([1.9,2.05])
% title(['$\xi = $',num2str(xi)],'Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$E$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\bar{\sigma}_t$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight