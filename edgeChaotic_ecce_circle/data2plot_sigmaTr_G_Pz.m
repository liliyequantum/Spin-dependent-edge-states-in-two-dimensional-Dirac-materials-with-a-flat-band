clear
close all
maxNumCompThreads(15)

%%%%%%%%%%% EE step size must consistence for Pz %%%%%%%%%%%%%%%%%%%%%

% edge, sigma_tr, G
load('./edge/data/data_sigmaDiffEdge_theta_E_xi.mat')
disp('input edge data')
oneMCosThetaEdge = 1 - cos(repmat(Theta,length(EE),1));

sigmaTrEdge = zeros(length(EE),length(Xi));
GEdge = zeros(length(EE),length(Xi));
for xiIdx = 1:1:length(Xi)
    sigmaDiff = sigmaDiffEdge(:,:,xiIdx);
    sigmaTrEdge(:,xiIdx) = sum(oneMCosThetaEdge.*sigmaDiff,2)*(Theta(2)-Theta(1));
    GEdge(:,xiIdx) = transpose(EE)./sigmaTrEdge(:,xiIdx);
    disp(xiIdx/length(Xi))
end

sigmaTrEdgeAveE = sum(sigmaTrEdge)*(EE(2)-EE(1))/(EE(end)-EE(1));
GEdgeAveE = sum(GEdge)*(EE(2)-EE(1))/(EE(end)-EE(1));

% chaotic, sigma_tr, G
mu = 20;
load(['./chaotic/data/data_sigmaDiffChaotic_theta_E_xi_mu',num2str(mu),'.mat'])
disp('input chaotic data')
oneMCosThetaChaotic = 1 - cos(repmat(Theta,length(EE),1));

sigmaTrChaotic = zeros(length(EE),length(Xi));
GChaotic = zeros(length(EE),length(Xi));
for xiIdx = 1:1:length(Xi)
    sigmaDiff = sigmaDiffChaotic(:,:,xiIdx);
    sigmaTrChaotic(:,xiIdx) = sum(oneMCosThetaChaotic.*sigmaDiff,2)*(Theta(2)-Theta(1));
    GChaotic(:,xiIdx) = transpose(EE)./sigmaTrChaotic(:,xiIdx);
    disp(xiIdx/length(Xi))
end

sigmaTrChaoticAveE = sum(sigmaTrChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));
GChaoticAveE = sum(GChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));

% Pz
pz = (sigmaTrEdge - sigmaTrChaotic)./(sigmaTrEdge + sigmaTrChaotic);
pzAveE = sum(pz)*(EE(2)-EE(1))/(EE(end)-EE(1));

pzG = (GChaotic - GEdge)./(GChaotic + GEdge);
pzGAveE = sum(pzG)*(EE(2)-EE(1))/(EE(end)-EE(1));

if exist('pictureData','dir')==0
    mkdir('pictureData');
end

save(['./pictureData/dataToPlot_sigmaTr_G_Pz_mu_',num2str(mu),'.mat'],'mu','Xi','sigmaTrEdgeAveE','sigmaTrChaoticAveE',...
    'GEdgeAveE','GChaoticAveE','pzAveE','pzGAveE')