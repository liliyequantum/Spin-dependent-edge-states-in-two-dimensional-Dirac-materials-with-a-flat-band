clear
close all
load('./pictureData/dataToPlot_sigmaTrOverTotal_Pz_xi_mu_homogen.mat')

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
pcolor(XiM, muM, pzOverTotalAveE); hold on;
title('$\langle q_{z} \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25,'Fontname', 'Times New Roman')
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight

f = figure();
plot(mmu,pzOverTotalAveE(:,1),'color',color_yellow,'LineWidth',3)
hold on
scatter(mmu,pzOverTotalAveE(:,1),10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',2.5)
hold on
% title('$\xi=0$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
set(gca, 'LineWidth',3,'Fontname', 'Times New Roman','FontSize',35)
xlabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',35)
ylabel('$\langle q_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',35)
axis tight
yticks([0.2,0.4,0.6,0.8,0.96])
xticks([20,24])
% ylim([0,1])
grid on

muNum = 11;
f = figure();
plot(Xi,pzOverTotalAveE(muNum,:),'color',color_yellow,'LineWidth',2)
hold on
scatter(Xi,pzOverTotalAveE(muNum,:),10,'MarkerEdgeColor', color_yellow,...
              'MarkerFaceColor', color_yellow,...
              'LineWidth',1.5)
hold on
title(['$\mu=$',num2str(mmu(muNum))],'Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\langle q_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
ylim([0,1])

maxPzOverTotalAveE = max(transpose(pzOverTotalAveE));
figure()
plot(mmu,maxPzOverTotalAveE,'color',color_purple,'LineWidth',2)
hold on
scatter(mmu,maxPzOverTotalAveE,10,'MarkerEdgeColor', color_purple,...
              'MarkerFaceColor', color_purple,...
              'LineWidth',1.5)
hold on
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\max\langle q_z \rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight
ylim([0,1])
% grid on

figure()
pcolor(XiM, muM, sigmaTrChaoticAveE); hold on;
title('$\langle \sigma^{\uparrow}_{tr}\rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25,'Fontname', 'Times New Roman')
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight

figure()
pcolor(XiM, muM, sigmaTotalChaoticAveE); hold on;
title('$\langle \sigma^{\uparrow}_{t}\rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25,'Fontname', 'Times New Roman')
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight

figure()
pcolor(XiM, muM, sigmaTrChaoticAveE./sigmaTotalChaoticAveE); hold on;
title('$\langle \sigma_{tr}\rangle/\langle \sigma_{t}\rangle$ Chaotic','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight

figure()
pcolor(XiM, muM, sigmaTrOverTotalChaoticAveE); hold on;
title('$\langle \sigma_{tr}/\sigma_{t}\rangle$ Chaotic','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
shading flat; 
colorbar;
set(gca,'fontsize',25)
xlabel('$\xi$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\mu$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight


