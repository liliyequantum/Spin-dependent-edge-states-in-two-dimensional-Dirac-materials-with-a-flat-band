clear
close all
maxNumCompThreads(15)

%%%%%%%%%%% EE step size must consistence for Pz %%%%%%%%%%%%%%%%%%%%%

% edge, sigma_tr
load('./edge/data/data_sigmaDiffEdge_theta_E_xi.mat')
disp('input edge data')
oneMCosThetaEdge = 1 - cos(repmat(Theta,length(EE),1));

sigmaTrEdge = zeros(length(EE),length(Xi));
sigmaTotalEdge = zeros(length(EE),length(Xi));
for xiIdx = 1:1:length(Xi)
    sigmaDiff = sigmaDiffEdge(:,:,xiIdx);
    sigmaTrEdge(:,xiIdx) = sum(oneMCosThetaEdge.*sigmaDiff,2)*(Theta(2)-Theta(1));
    sigmaTotalEdge(:,xiIdx) = sum(sigmaDiff,2)*(Theta(2)-Theta(1));
    disp(xiIdx/length(Xi))
end

sigmaTrTotalEdgeAveE = sum(sigmaTrEdge./sigmaTotalEdge)*(EE(2)-EE(1))/(EE(end)-EE(1));

% chaotic, sigma_tr, G
mu = 24;
load(['./chaotic/data/data_sigmaDiffChaotic_theta_E_xi_mu',num2str(mu),'.mat'])
disp('input chaotic data')
oneMCosThetaChaotic = 1 - cos(repmat(Theta,length(EE),1));

sigmaTrChaotic = zeros(length(EE),length(Xi));
sigmaTotalChaotic = zeros(length(EE),length(Xi));
for xiIdx = 1:1:length(Xi)
    sigmaDiff = sigmaDiffChaotic(:,:,xiIdx);
    sigmaTrChaotic(:,xiIdx) = sum(oneMCosThetaChaotic.*sigmaDiff,2)*(Theta(2)-Theta(1));
    sigmaTotalChaotic(:,xiIdx) = sum(sigmaDiff,2)*(Theta(2)-Theta(1));
    disp(xiIdx/length(Xi))
end

sigmaTrTotalChaoticAveE = sum(sigmaTrChaotic./sigmaTotalChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));

% Pz
pzOverTotal = (sigmaTrEdge./sigmaTotalEdge - sigmaTrChaotic./sigmaTotalChaotic)...
    ./(sigmaTrEdge./sigmaTotalEdge + sigmaTrChaotic./sigmaTotalChaotic);
pzOverTotalAveE = sum(pzOverTotal)*(EE(2)-EE(1))/(EE(end)-EE(1));


if exist('pictureData','dir')==0
    mkdir('pictureData');
end

save(['./pictureData/dataToPlot_sigmaTrOverTotal_Pz_mu_',num2str(mu),'.mat'],...
    'mu','Xi','sigmaTrTotalEdgeAveE','sigmaTrTotalChaoticAveE',...
    'pzOverTotalAveE')