clear
close all

dbstop if warning
% dbclear if warning

xi = 0;
spinIndex = 1; % up = 1, down = -1
incidentAngle = 0;

if spinIndex==-1
    v1 = -10; 
    v2 = 40;
    L = -40:1:40;
    EE = [14.5805,14.8,16];
else
    mu = 24;
    v1 = -10 - 2*mu;
    v2 = 40 - 2*mu;
    L = -35:1:35;
    EE = [14.8,15.075,16];
    % EE = [14.8,15.0235,16];
end

R_1 = 1;
R_2 = 0.6;

Theta = linspace(-pi,pi,1000); %scattering theta
sigmaDiff = zeros(length(EE),length(Theta));
oneMinusCosSigmaDiff = zeros(length(EE),length(Theta));

A_0 = diag(0.5/sqrt(2)*1j.^(L-1).*exp(-1i*incidentAngle*L));
a0row = transpose(diag(A_0)); 
for NN=1:1:length(EE)
    E = EE(NN);
    
    k_0 = abs(E);
    k_1 = abs(E - v1);
    k_2 = abs(E - v2);
    alpha_0 = sign(E);
    alpha_1 = sign(E - v1);
    alpha_2 = sign(E - v2);
    
    X_1 = diag(besselh(L, 1, k_0*R_1));
    X_2 = diag(besselh(L, 2, k_0*R_1));
    x_1 = diag(besselh(L, 1, k_1*R_1));
    x_2 = diag(besselh(L, 2, k_1*R_1));
    
    Y_1 = diag(besselh(L+1, 1, k_0*R_1));
    Y_2 = diag(besselh(L+1, 2, k_0*R_1));
    y_1 = diag(besselh(L+1, 1, k_1*R_1));
    y_2 = diag(besselh(L+1, 2, k_1*R_1));
    
    Z_1 = diag(besselh(L-1, 1, k_0*R_1));
    Z_2 = diag(besselh(L-1, 2, k_0*R_1));
    z_1 = diag(besselh(L-1, 1, k_1*R_1));
    z_2 = diag(besselh(L-1, 2, k_1*R_1));
    
    temp_1 = besselj(L, k_2*R_2).*(besselh(L-1, 2, k_1*R_2) -besselh(L+1, 2, k_1*R_2)) ...
             - alpha_1*alpha_2*besselh(L, 2, k_1*R_2).*(besselj(L-1, k_2*R_2) - besselj(L+1, k_2*R_2));
    temp_2 =   besselj(L, k_2*R_2).*(besselh(L-1, 1, k_1*R_2) -besselh(L+1, 1, k_1*R_2)) ...
             - alpha_1*alpha_2*besselh(L, 1, k_1*R_2).*(besselj(L-1, k_2*R_2) - besselj(L+1, k_2*R_2));
    S_cd = diag(-temp_1./temp_2);
    [column,row] = meshgrid(L,L);
    
    
    U = besselj(column - row, k_1*xi);
    U_inv = besselj(row - column, k_1*xi);
    S_od = U_inv*S_cd*U;
    
    F = x_2 + S_od*x_1;
    G = y_2 + S_od*y_1;
    H = z_2 + S_od*z_1;
    T = F\(H-G);
    
    S = -(Z_2 - Y_2-alpha_0*alpha_1*X_2*T)/(Z_1 - Y_1 - alpha_0*alpha_1*X_1*T);
    TM = S - eye(size(S));
    
    for thetaIdx = 1:1:length(Theta)
        theta = Theta(thetaIdx);
        thetacol = transpose((-1i).^L.*exp(1i*L*theta));
        oneMinusCosSigmaDiff(NN,thetaIdx) = (1-cos(theta))*abs(a0row*TM*thetacol)^2*2/pi/k_0;
        sigmaDiff(NN,thetaIdx) = abs(a0row*TM*thetacol)^2*2/pi/k_0;
    end
   
    disp(NN/length(EE))
end       

%% plot

color_blue = [31, 119, 180]/255;
color_red = [214, 39, 40]/255;
color_yellow = [255, 127, 14]/255;
color_purple = [148, 103, 189]/255;
color_green = [44, 160, 44]/255;

if length(EE)==1
    red2yellow = color_red;
    blue2green = color_blue;
else
    numcolor = length(EE);
    red2yellow = [0, 0, 0;255,64,64;247,250,19]./255;%colorGradient([233,66,53]/255,[250,187,5]/255,numcolor);
    blue2green = [64,64,255;0,0,0;42,245,152]./255;%colorGradient(color_blue,[42,245,152]/255,numcolor);
end

for NN = 1:1:length(EE)  
    f = figure(1);
    f.Position = [100 100 800 400];
    if spinIndex==1
        plot(Theta/pi,sigmaDiff(NN,:),'color',red2yellow(NN,:),'LineWidth',2)
        hold on
    else
        plot(Theta/pi,sigmaDiff(NN,:),'color',blue2green(NN,:),'LineWidth',2)
        hold on
    end    
end
%     ylim([0,2])
xlabel('$\theta$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',20)
xlabel('$\theta(\pi)$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
ylabel('$\sigma_{diff}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',20)
axis tight
temp = cell(1,length(EE));
for i = 1:1:length(EE)
    temp{1,i} = ['$E = $',num2str(EE(i))];
end
legend(temp,'interpreter','latex','fontsize',20)

for NN = 1:1:length(EE) 
    f = figure(2);
    f.Position = [100 100 800 400];
    if spinIndex==1
        plot(Theta/pi,oneMinusCosSigmaDiff(NN,:),'color',red2yellow(NN,:),'LineWidth',3)
        hold on
    else
        plot(Theta/pi,oneMinusCosSigmaDiff(NN,:),'color',blue2green(NN,:),'LineWidth',3)
        hold on
    end    
end
temp = cell(1,length(EE));
for i = 1:1:length(EE)
    if i==2
        temp{1,i} = [''];
    else
        temp{1,i} = ['$E = $',num2str(EE(i))];
    end
end
legend(temp,'interpreter','latex','fontsize',25)
% yticks([0.5,1])
yticks([0.05,0.1])
%     ylim([0,2])
set(gca, 'LineWidth',1,'Fontname', 'Times New Roman','FontSize',25)
xlabel('$\theta(\pi)$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
ylabel('$f(\theta)\sigma_{\rm diff}$','Interpreter','latex','Fontname', 'Times New Roman','FontSize',25)
axis tight

% ylim([0,0.02])

