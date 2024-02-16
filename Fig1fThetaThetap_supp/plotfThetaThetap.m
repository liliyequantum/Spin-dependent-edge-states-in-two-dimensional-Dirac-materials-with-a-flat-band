clear
close all

load('fthetathetap_mu_22_supp.mat')
%% plot

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

if length(Xi)==1
    red2yellow = color_red;
    blue2green = color_blue;
else
    numcolor = length(Xi);
    red2yellow = [[42,245,152]/255;[240,5,5]/255;[254,240,1]/255];%colorGradient([240,5,5]/255,[254,240,1]/255,numcolor);
    blue2green = colorGradient([8,179,229]/255,[42,245,152]/255,numcolor);
end

% for xiIdx = 1:1:length(Xi)  
%     f = figure(1);
%     f.Position = [100 100 600 400];
%     if spinIndex==1
%         plot(Theta,sigmaDiffAveE(1,:,xiIdx),'color',red2yellow(xiIdx,:),'LineWidth',2)
%         hold on
%     else
%         plot(Theta,sigmaDiffAveE(1,:,xiIdx),'color',blue2green(xiIdx,:),'LineWidth',2)
%         hold on
%     end    
% end
% %     ylim([0,2])
% xlabel('$\theta$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% title(['v1 = ',num2str(v1),' v2 = ',num2str(v2)],'interpreter','latex','fontsize',20)
% set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
% xlabel('$\theta$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% ylabel('$\langle\sigma_{diff}\rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
% axis tight
% temp = cell(1,length(Xi));
% for i = 1:1:length(Xi)
%     temp{1,i} = ['$\xi = $',num2str(Xi(i))];
% end
% legend(temp,'interpreter','latex','fontsize',20)

for xiIdx = 1:1:length(Xi)  
    f = figure(2);
    f.Position = [100 100 800 400];
    if spinIndex==1
        plot(Theta/pi,oneMinusCosSigmaDiffAveE(1,:,xiIdx),'color',red2yellow(xiIdx,:),'LineWidth',3)
        hold on
    else
        plot(Theta/pi,oneMinusCosSigmaDiffAveE(1,:,xiIdx),'color',blue2green(xiIdx,:),'LineWidth',3)
        hold on
    end    
end
%     ylim([0,2])
% title(['v1 = ',num2str(v1),' v2 = ',num2str(v2)],'interpreter','latex','fontsize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',30)
xlabel('$\theta(\pi)$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
ylabel('$f(\theta)\langle\sigma_{diff}\rangle$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',30)
axis tight
temp = cell(1,length(Xi));
for i = 1:1:length(Xi)
    temp{1,i} = ['$\xi = $',num2str(Xi(i))];
end
legend(temp,'interpreter','latex','fontsize',25)
set(gca,'box','off')
yticks([0.1,0.2])