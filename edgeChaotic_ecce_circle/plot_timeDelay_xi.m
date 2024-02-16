clear
close all

load('./chaotic/data/data_timeDelayUp_xi_mu_24.mat')
load('./edge/data/data_timeDelayDown_xi.mat')

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

time_delayUpMax = max(real(tauMatrixUp));
time_delayDownMax = max(real(tauMatrixDown));

f = figure(1);
f.Position = [100 100 800 400];
semilogy(Xi,time_delayUpMax,'color',color_red,'LineWidth',2)
hold on
semilogy(Xi,time_delayDownMax,'color',color_blue,'LineWidth',2)
hold on
scatter(Xi,time_delayUpMax,10,'MarkerEdgeColor',[214, 39, 40]/255,...
              'MarkerFaceColor',[214, 39, 40]/255,...
              'LineWidth',1.5)
hold on
scatter(Xi,time_delayDownMax,10,'MarkerEdgeColor',[31, 119, 180]/255,...
              'MarkerFaceColor',[31, 119, 180]/255,...
              'LineWidth',1.5)
hold on
yticks([1e2,1e4])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
ylabel('$\tau_{max}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
grid on
% %%
% clear
% close all
% 
% color_blue = [31, 119, 180]/255;
% color_red = [214, 39, 40]/255;
% color_yellow = [255, 127, 14]/255;
% color_purple = [148, 103, 189]/255;
% color_green = [44, 160, 44]/255;
% 
% load('./chaotic/data/data_timeDelayUp_xi.mat')
% rtauMatrixUp = real(tauMatrixUp);
% 
% for i = 1:10:length(Xi)
%     figure(i)
%     plot(EE(1:end-1),rtauMatrixUp(:,i),'color',color_red,'LineWidth',2)
%     hold on
%     title(['spin up $\xi$ = ',num2str(Xi(i))],'interpreter','latex')
%     set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
%     xlabel('$E$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
%     ylabel('$\tau$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
%     axis tight
% end
% 
% [XiM,EM] = meshgrid(Xi,EE(1:end-1));
% figure()
% pcolor(XiM, EM, rtauMatrixUp); hold on;
% colormap(jet(20));
% shading flat; 
% colorbar;
% set(gca,'fontsize',15)

%%
% clear
% close all
% 
% color_blue = [31, 119, 180]/255;
% color_red = [214, 39, 40]/255;
% color_yellow = [255, 127, 14]/255;
% color_purple = [148, 103, 189]/255;
% color_green = [44, 160, 44]/255;
% 
% load('./edge/data/data_timeDelayDown_xi.mat')
% rtauMatrixDown = real(tauMatrixDown);
% % for i = 1:10:length(Xi)
% %     figure(length(Xi)+i)
% %     plot(EE(1:end-1),rtauMatrixDown(:,i),'color',color_blue,'LineWidth',2)
% %     hold on
% %     title(['spin up $\xi$ = ',num2str(Xi(i))],'interpreter','latex')
% %     set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
% %     xlabel('$E$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% %     ylabel('$\tau$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% %     axis tight
% % end
% 
% [XiM,EM] = meshgrid(Xi,EE(1:end-1));
% figure()
% pcolor(XiM, EM, rtauMatrixDown); hold on;
% % colormap(jet(20));
% shading flat; 
% colorbar;
% set(gca,'fontsize',15)