function RelationwavefunctionAndCrossSection(mu,xi,E,N,xleft,xright)
    maxNumCompThreads(15)
    
%     mu = 24;
%     xi = 0;
%     E = 14.8;
% 
%     N=500;
%     xleft = 1;%1.1;
%     xright = 6;

    coeForWave(mu,xi,E)
    
    load(['./data/coeDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E)])

    % total region

    x_choose=linspace(xleft,xright,N);
    y_choose=linspace(-3,3,N);
    dx = x_choose(2)-x_choose(1);
    dy = y_choose(2)-y_choose(1);

    [xx,yy]=meshgrid(x_choose,y_choose);
    zz=xx+sqrt(-1)*yy;

    r = abs(zz);
    theta = angle(zz);
    rp = abs(zz-xi);
    thetap = angle(zz-xi);
    
    psiA = zeros(N,N);
    psiB = zeros(N,N);
    psiC = zeros(N,N);
    % psi_I
    % disp('psi_I, o origin')
    TM = S - eye(size(S));
    for mm = 1:1:length(L)
            for ll = 1:1:length(L)                
                psiA = psiA + A_0(mm,mm)*TM(mm,ll)*besselh(L(ll)-1, 1, k_0*r).*exp(1j*(L(ll)-1)*theta);
                psiB = psiB + A_0(mm,mm)*TM(mm,ll)*1j*alpha_0*sqrt(2)*besselh(L(ll), 1, k_0*r).*exp(1j*L(ll)*theta);
                psiC = psiC - A_0(mm,mm)*TM(mm,ll)*besselh(L(ll)+1, 1, k_0*r).*exp(1j*(L(ll)+1)*theta);
            end
           disp({'scattering I, o: ', num2str(mm/length(L))})
    end
 
    psiA = psiA/sqrt(2);
    psiB = psiB/sqrt(2);
    psiC = psiC/sqrt(2);
    
    probDens = abs(psiA).^2+abs(psiB).^2+abs(psiC).^2;
    save(['./data/data2plotDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_N_',num2str(N),'_xleft_',...
        num2str(xleft),'_xright_',num2str(xright),'.mat'],...
        'xx','yy','probDens')
    
end
