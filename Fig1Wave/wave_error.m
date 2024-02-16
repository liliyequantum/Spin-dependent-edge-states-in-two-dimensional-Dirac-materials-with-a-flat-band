function wave_error(mu,xi,E,Lrange,N,xyra)
    maxNumCompThreads(15)
    
    coeForWave(mu,xi,E,Lrange)
    
    load(['./data/coeDensity_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'.mat'])
    
    % boundary1
    
    Nj1=800;
    theta1=linspace(0,2*pi,Nj1+1);
    theta1=theta1(1:Nj1);
    Zj1=R_1*exp(sqrt(-1)*theta1);
    
    % boundary2
    
    Nj2=800;
    theta2op=linspace(0,2*pi,Nj2+1);
    theta2op=theta2op(1:Nj2);
    Zj2op=R_2*exp(sqrt(-1)*theta2op);
    Zj2 = xi + Zj2op;
    
    % total region
    
%     N=450;
%     xyra = 1.1;

    x_choose=linspace(-xyra,xyra,N);
    y_choose=linspace(-xyra,xyra,N);
    dx = x_choose(2)-x_choose(1);
    dy = y_choose(2)-y_choose(1);

    [xx,yy]=meshgrid(x_choose,y_choose);
    zz=xx+sqrt(-1)*yy;
    
    In1=zeros(N,N);
    In2=zeros(N,N);
    In3=zeros(N,N);
    for i=1:N
        for j=1:N
            if inpolygon(xx(i,j),yy(i,j),real(Zj2),imag(Zj2))
                In3(i,j)=1;
            elseif inpolygon(xx(i,j),yy(i,j),real(Zj1),imag(Zj1))
                In2(i,j)=1;
            else
                In1(i,j)=1;
            end
        end
        disp(i/N)
    end
    
    r = abs(zz);
    theta = angle(zz);
    rp = abs(zz-xi);
    thetap = angle(zz-xi);
    
    psiA = zeros(N,N);
    psiB = zeros(N,N);
    psiC = zeros(N,N);
    % psi_I
    % disp('psi_I, o origin')
    for mm = 1:1:length(L)
            for ll = 1:1:length(L)
                if ll==mm
                    psiA = psiA + A_0(mm,mm)*besselh(L(mm)-1, 2, k_0*r).*exp(1j*(L(mm)-1)*theta).*In1;
                    psiB = psiB + A_0(mm,mm)*1j*alpha_0*sqrt(2)*besselh(L(mm), 2, k_0*r).*exp(1j*L(mm)*theta).*In1;
                    psiC = psiC - A_0(mm,mm)*besselh(L(mm)+1, 2, k_0*r).*exp(1j*(L(mm)+1)*theta).*In1;
                end
                psiA = psiA + A_0(mm,mm)*S(mm,ll)*besselh(L(ll)-1, 1, k_0*r).*exp(1j*(L(ll)-1)*theta).*In1;
                psiB = psiB + A_0(mm,mm)*S(mm,ll)*1j*alpha_0*sqrt(2)*besselh(L(ll), 1, k_0*r).*exp(1j*L(ll)*theta).*In1;
                psiC = psiC - A_0(mm,mm)*S(mm,ll)*besselh(L(ll)+1, 1, k_0*r).*exp(1j*(L(ll)+1)*theta).*In1;
            end
           disp({'psi_I, o: ', num2str(mm/length(L))})
    end
    
    % psi_II 
    % disp('psi_II, oprime origin')
    for mm = 1:1:length(L)
        for ll = 1:1:length(L)
    
            psiA = psiA + A_tilde(mm,ll)*(besselh(L(ll)-1, 2, k_1*rp).*exp(1j*(L(ll)-1)*thetap)+...
                S_cd(ll,ll)*besselh(L(ll)-1, 1, k_1*rp).*exp(1j*(L(ll)-1)*thetap)).*In2;
            psiB = psiB + A_tilde(mm,ll)*1j*alpha_1*sqrt(2)*(besselh(L(ll), 2, k_1*rp).*exp(1j*L(ll)*thetap)+...
                S_cd(ll,ll)*besselh(L(ll), 1, k_1*rp).*exp(1j*L(ll)*thetap)).*In2;
            psiC = psiC - A_tilde(mm,ll)*(besselh(L(ll)+1, 2, k_1*rp).*exp(1j*(L(ll)+1)*thetap)+...
                S_cd(ll,ll)*besselh(L(ll)+1, 1, k_1*rp).*exp(1j*(L(ll)+1)*thetap)).*In2;
    
        end
        disp({'psi_II, oprime: ', num2str(mm/length(L))})
    end
    
    
    % psi_III
    % disp('psi_III, oprime origin')
    for mm = 1:1:length(L)
        for ll = 1:1:length(L)
            psiA = psiA + B_tilde(mm,ll)*besselj(L(ll)-1, k_2*rp).*exp(1j*(L(ll)-1)*thetap).*In3;
            psiB = psiB + B_tilde(mm,ll)*1j*alpha_2*sqrt(2)*besselj(L(ll), k_2*rp).*exp(1j*L(ll)*thetap).*In3;
            psiC = psiC - B_tilde(mm,ll)*besselj(L(ll)+1, k_2*rp).*exp(1j*(L(ll)+1)*thetap).*In3;
        end
        disp({'psi_III, oprime: ', num2str(mm/length(L))})
    end
    
    psiA = psiA/sqrt(2);
    psiB = psiB/sqrt(2);
    psiC = psiC/sqrt(2);
    
    probDens = abs(psiA).^2+abs(psiB).^2+abs(psiC).^2;    
    save(['./data/data2plotDensityAndError_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'],...
        'xx','yy','probDens','Zj1','Zj2','In1','In2','In3')
   
    %% error
    disp('equation error')
    
    dx=xx(1,2)-xx(1,1);
    dy=yy(2,1)-yy(1,1);
    err=zeros(N,N);
    
    
    for i=2:N-1
        for j=2:N-1
            
            partial_Ax=(psiA(i,j+1)-psiA(i,j-1))/(2*dx);
            partial_Bx=(psiB(i,j+1)-psiB(i,j-1))/(2*dx);
            partial_Cx=(psiC(i,j+1)-psiC(i,j-1))/(2*dx);
            partial_Ay=(psiA(i+1,j)-psiA(i-1,j))/(2*dy);
            partial_By=(psiB(i+1,j)-psiB(i-1,j))/(2*dy);
            partial_Cy=(psiC(i+1,j)-psiC(i-1,j))/(2*dy);
            
            if In1(i-1,j)==1 && In1(i+1,j)==1 && In1(i,j+1)==1 && In1(i,j-1)==1
                
                left1=-sqrt(-1)*partial_Bx-partial_By;
                left2=-sqrt(-1)*partial_Ax-sqrt(-1)*partial_Cx+partial_Ay-partial_Cy;
                left3=-sqrt(-1)*partial_Bx+partial_By;
                right1=sqrt(2)*k_0*alpha_0*psiA(i,j);
                right2=sqrt(2)*k_0*alpha_0*psiB(i,j);
                right3=sqrt(2)*k_0*alpha_0*psiC(i,j);
                err(i,j)=(abs(left1-right1)+abs(left2-right2)+abs(left3-right3))/...
                    (abs(left1)+abs(left2)+abs(left3)+abs(right1)+abs(right2)+abs(right3)); % definition?
                
            elseif In2(i-1,j)==1 && In2(i+1,j)==1 && In2(i,j+1)==1 && In2(i,j-1)==1
                
                left1=-sqrt(-1)*partial_Bx-partial_By;
                left2=-sqrt(-1)*partial_Ax-sqrt(-1)*partial_Cx+partial_Ay-partial_Cy;
                left3=-sqrt(-1)*partial_Bx+partial_By;
                right1=sqrt(2)*k_1*alpha_1*psiA(i,j);
                right2=sqrt(2)*k_1*alpha_1*psiB(i,j);
                right3=sqrt(2)*k_1*alpha_1*psiC(i,j);
                err(i,j)=(abs(left1-right1)+abs(left2-right2)+abs(left3-right3))/...
                    (abs(left1)+abs(left2)+abs(left3)+abs(right1)+abs(right2)+abs(right3)); % definition?
                
            elseif In3(i-1,j)==1 && In3(i+1,j)==1 && In3(i,j+1)==1 && In3(i,j-1)==1
                
                left1=-sqrt(-1)*partial_Bx-partial_By;
                left2=-sqrt(-1)*partial_Ax-sqrt(-1)*partial_Cx+partial_Ay-partial_Cy;
                left3=-sqrt(-1)*partial_Bx+partial_By;
                right1=sqrt(2)*k_2*alpha_2*psiA(i,j);
                right2=sqrt(2)*k_2*alpha_2*psiB(i,j);
                right3=sqrt(2)*k_2*alpha_2*psiC(i,j);
                err(i,j)=(abs(left1-right1)+abs(left2-right2)+abs(left3-right3))/...
                    (abs(left1)+abs(left2)+abs(left3)+abs(right1)+abs(right2)+abs(right3)); % definition?
                
            end
        end
        disp(['progress ', num2str(i/N)])
    end
    
    save(['./data/data2plotDensityAndError_mu_',num2str(mu),'_xi_',...
        num2str(xi),'_E_',num2str(E),'_L_',num2str(Lrange),'_N_',num2str(N),'_xyra_',...
        num2str(xyra),'.mat'],...
        'xx','yy','probDens','err','Zj1','Zj2','In1','In2','In3')
    
end
