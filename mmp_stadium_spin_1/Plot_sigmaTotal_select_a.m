clear
close all

load('data2Plot_sigmaTr_E_selecta.mat')

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

aaEdge = 0:0.05:1.5;
aaChaotic = 0:0.01:1.5;
aIdxSelectChaotic = [1,51,101]; % a= 0.5,1

blue2green = colorGradient([8,179,229]/255,[42,245,152]/255,length(aIdxSelectEdge));
red2yellow = colorGradient([240,5,5]/255,[254,240,1]/255,length(aIdxSelectChaotic));

f = figure(1);
f.Position = [100 100 600 400];
% f.Position = [100 100 600 400];
for aIdx = 1:1:length(aIdxSelectEdge)
    plot(EEEdge*R, transpose(sigmaTotalEdge(:,aIdx)),'color',blue2green(aIdx,:),...
        'LineWidth',2,'DisplayName',['a = ',num2str(aaEdge(aIdxSelectEdge(aIdx)))])
    hold on
end
legend('FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$ER$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\sigma_{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylim([90,120])

f = figure(2);
f.Position = [100 100 600 400];
for aIdx = 1:1:length(aIdxSelectChaotic)
    plot(EEChaotic*R, transpose(sigmaTotalChaotic(:,aIdx)),'color',red2yellow(aIdx,:),...
        'LineWidth',2,'DisplayName',['a = ',num2str(aaChaotic(aIdxSelectChaotic(aIdx)))])
    hold on
end
legend('FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$ER$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$\sigma_{t}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
