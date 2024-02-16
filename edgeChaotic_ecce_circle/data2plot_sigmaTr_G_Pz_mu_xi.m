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
end

disp('To input chaotic data')
% chaotic, sigma_tr, G
a = 17.5:0.25:23;
c = 24.75:0.25:27.5;
b = 23+0.0625:0.0625:24.75-0.0625;
mmu = [a,b,c];
pzAveE = zeros(length(mmu),length(Xi));
pzGAveE = zeros(length(mmu),length(Xi));
sigmaTrChaoticAveE = zeros(length(mmu),length(Xi));
GChaoticAveE = zeros(length(mmu),length(Xi));
for muIdx = 1:1:length(mmu)
    mu = mmu(muIdx);
    load(['./chaotic/data/data_sigmaDiffChaotic_theta_E_xi_mu',num2str(mu),'.mat'])
    oneMCosThetaChaotic = 1 - cos(repmat(Theta,length(EE),1));
    
    sigmaTrChaotic = zeros(length(EE),length(Xi));
    GChaotic = zeros(length(EE),length(Xi));
    for xiIdx = 1:1:length(Xi)
        sigmaDiff = sigmaDiffChaotic(:,:,xiIdx);
        sigmaTrChaotic(:,xiIdx) = sum(oneMCosThetaChaotic.*sigmaDiff,2)*(Theta(2)-Theta(1));
        GChaotic(:,xiIdx) = transpose(EE)./sigmaTrChaotic(:,xiIdx);
    end

    sigmaTrChaoticAveE(muIdx,:) = sum(sigmaTrChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));
    GChaoticAveE(muIdx,:) = sum(GChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));
    
    % Pz
    pz = (sigmaTrEdge - sigmaTrChaotic)./(sigmaTrEdge + sigmaTrChaotic);
    pzAveE(muIdx,:) = sum(pz)*(EE(2)-EE(1))/(EE(end)-EE(1));

    pzG = (GChaotic - GEdge)./(GChaotic + GEdge);
    pzGAveE(muIdx,:) = sum(pzG)*(EE(2)-EE(1))/(EE(end)-EE(1));

    disp(muIdx/length(mmu))

end

if exist('pictureData','dir')==0
    mkdir('pictureData');
end

save('./pictureData/dataToPlot_sigmaTr_G_Pz_mu17p5to27p5MoreFine_xi.mat','Xi','mmu',"pzAveE",'pzGAveE',...
    'GChaoticAveE','sigmaTrChaoticAveE')

