function inOutwave(mu,xi,E,Lrange,N,xyra)
    maxNumCompThreads(15)
    
    coeForWave(mu,xi,E,Lrange)
    
    load(['./data/coeDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'.mat'])

    % region I boundary
    Rj=1;
    Nj=800;
    theta=linspace(0,2*pi,Nj+1);
    theta=theta(1:Nj);
    Zj=Rj*exp(sqrt(-1)*theta);
    
%     N=600;
%     xyra = 1.5;
    x_choose=linspace(-xyra,xyra,N);
    y_choose=linspace(-xyra,xyra,N);
    dx = x_choose(2)-x_choose(1);
    dy = y_choose(2)-y_choose(1);

    [xx,yy]=meshgrid(x_choose,y_choose);
    zz=xx+sqrt(-1)*yy;
    
    Out=ones(N,N);
    for i=1:N
        for j=1:N
            if inpolygon(xx(i,j),yy(i,j),real(Zj),imag(Zj))   
                Out(i,j)=0;
            end
        end
        disp(i/N)
    end

    r = abs(zz);
    theta = angle(zz);
    rp = abs(zz-xi);
    thetap = angle(zz-xi);
    
    psiAIn = zeros(N,N);
    psiBIn = zeros(N,N);
    psiCIn = zeros(N,N);
   
    % psiIn
    for mm = 1:1:length(L)
        psiAIn = psiAIn + 2*A_0(mm,mm)*besselj(L(mm)-1, k_0*r).*exp(1j*(L(mm)-1)*theta);
        psiBIn = psiBIn + 2*A_0(mm,mm)*1j*alpha_0*sqrt(2)*besselj(L(mm), k_0*r).*exp(1j*L(mm)*theta);
        psiCIn = psiCIn - 2*A_0(mm,mm)*besselj(L(mm)+1, k_0*r).*exp(1j*(L(mm)+1)*theta);
         disp({'psiIn: ', num2str(mm/length(L))})
    end

    psiAIn = psiAIn/sqrt(2).*Out;
    psiBIn = psiBIn/sqrt(2).*Out;
    psiCIn = psiCIn/sqrt(2).*Out;
    probIn = abs(psiAIn).^2+abs(psiBIn).^2+abs(psiCIn).^2;

    %psiOut
    psiAOut = zeros(N,N);
    psiBOut = zeros(N,N);
    psiCOut = zeros(N,N);
    TM = S - eye(size(S));
    for mm = 1:1:length(L)
            for ll = 1:1:length(L)
                psiAOut = psiAOut + A_0(mm,mm)*TM(mm,ll)*besselh(L(ll)-1, 1, k_0*r).*exp(1j*(L(ll)-1)*theta);
                psiBOut = psiBOut + A_0(mm,mm)*TM(mm,ll)*1j*alpha_0*sqrt(2)*besselh(L(ll), 1, k_0*r).*exp(1j*L(ll)*theta);
                psiCOut = psiCOut - A_0(mm,mm)*TM(mm,ll)*besselh(L(ll)+1, 1, k_0*r).*exp(1j*(L(ll)+1)*theta);
            end
           disp({'psiOut: ', num2str(mm/length(L))})
    end
    psiAOut = psiAOut/sqrt(2).*Out;
    psiBOut = psiBOut/sqrt(2).*Out;
    psiCOut = psiCOut/sqrt(2).*Out;
    probOut = abs(psiAOut).^2+abs(psiBOut).^2+abs(psiCOut).^2;
   
    psiA = psiAIn + psiAOut;
    psiB = psiBIn + psiBOut;
    psiC = psiCIn + psiCOut;
    prob = abs(psiA).^2+abs(psiB).^2+abs(psiC).^2;

    save(['./data/data2plotInOutWave_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'],...
        'xx','yy','Out','Zj','probIn','probOut','prob')
    
end
