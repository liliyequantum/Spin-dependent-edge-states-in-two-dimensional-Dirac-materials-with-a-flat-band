clear
close all
maxNumCompThreads(15)
%%%%%%%%%%% EE step size must consistence for Pz %%%%%%%%%%%%%%%%%%%%%

% edge, sigma_tr, G
load('./edge/data/data_sigmaDiffEdge_theta_E_xi.mat')
disp('input edge data')
oneMCosThetaEdge = 1 - cos(repmat(Theta,length(EE),1));

sigmaTrEdge = zeros(length(EE),length(Xi));
sigmaTotalEdge = zeros(length(EE),length(Xi));
for xiIdx = 1:1:length(Xi)
    sigmaDiff = sigmaDiffEdge(:,:,xiIdx);
    sigmaTrEdge(:,xiIdx) = sum(oneMCosThetaEdge.*sigmaDiff,2)*(Theta(2)-Theta(1));
    sigmaTotalEdge(:,xiIdx) = sum(sigmaDiff,2)*(Theta(2)-Theta(1));
end

disp('To input chaotic data')
% chaotic, sigma_tr, G
% a = 17.5:0.25:23;
% c = 24.75:0.25:27.5;
% b = 23+0.0625:0.0625:24.75-0.0625;
mmu = 17.5:0.25:27.5;
pzOverTotalAveE = zeros(length(mmu),length(Xi));
sigmaTrChaoticAveE = zeros(length(mmu),length(Xi));
sigmaTotalChaoticAveE = zeros(length(mmu),length(Xi));
sigmaTrOverTotalChaoticAveE = zeros(length(mmu),length(Xi));
for muIdx = 1:1:length(mmu)
    mu = mmu(muIdx);
    load(['./chaotic/data/data_sigmaDiffChaotic_theta_E_xi_mu',num2str(mu),'.mat'])
    oneMCosThetaChaotic = 1 - cos(repmat(Theta,length(EE),1));
    
    sigmaTrChaotic = zeros(length(EE),length(Xi));
    sigmaTotalChaotic = zeros(length(EE),length(Xi));
    sigmaTrOverTotalChaotic = zeros(length(EE),length(Xi));
    for xiIdx = 1:1:length(Xi)
        sigmaDiff = sigmaDiffChaotic(:,:,xiIdx);
        sigmaTrChaotic(:,xiIdx) = sum(oneMCosThetaChaotic.*sigmaDiff,2)*(Theta(2)-Theta(1));
        sigmaTotalChaotic(:,xiIdx) = sum(sigmaDiff,2)*(Theta(2)-Theta(1));
        sigmaTrOverTotalChaotic(:,xiIdx) = sigmaTrChaotic(:,xiIdx)./sigmaTotalChaotic(:,xiIdx);
    end

    sigmaTrChaoticAveE(muIdx,:) = sum(sigmaTrChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));
    sigmaTotalChaoticAveE(muIdx,:) = sum(sigmaTotalChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));
    sigmaTrOverTotalChaoticAveE(muIdx,:) = sum(sigmaTrOverTotalChaotic)*(EE(2)-EE(1))/(EE(end)-EE(1));
    
    % Pz
    pz = (sigmaTrEdge./sigmaTotalEdge - sigmaTrChaotic./sigmaTotalChaotic)...
        ./(sigmaTrEdge./sigmaTotalEdge + sigmaTrChaotic./sigmaTotalChaotic);
    pzOverTotalAveE(muIdx,:) = sum(pz)*(EE(2)-EE(1))/(EE(end)-EE(1));

    disp(muIdx/length(mmu))

end

if exist('pictureData','dir')==0
    mkdir('pictureData');
end

save('./pictureData/dataToPlot_sigmaTrOverTotal_Pz_xi_mu_homogen.mat',...
    'Xi','mmu',"pzOverTotalAveE",...
    'sigmaTrChaoticAveE','sigmaTotalChaoticAveE','sigmaTrOverTotalChaoticAveE')

