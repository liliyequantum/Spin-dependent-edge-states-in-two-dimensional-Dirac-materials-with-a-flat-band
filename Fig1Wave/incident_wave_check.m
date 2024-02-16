clear
close all

mu = 0;
v1 = -10 - 2*mu;
v2 = 40 - 2*mu;
L = -40:1:40;
xi = 0;
E = 14.8;

R_1 = 1;
R_2 = 0.6;

A_0 = diag(0.5/sqrt(2)*1j.^(L-1));

%%
Rj=1;
Nj=300;
theta=linspace(0,2*pi,Nj+1);
theta=theta(1:Nj);
Zj=Rj*exp(sqrt(-1)*theta);

N=200; % 2000;
xyra = 1.5;
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

k_0 = abs(E);
alpha_0 = sign(E);
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
ps = real(psiAIn).^2+imag(psiAIn).^2+real(psiBIn).^2+imag(psiBIn).^2+real(psiCIn).^2+imag(psiCIn).^2;

figure()
pcolor(xx, yy, real(psiAIn)); hold on;
title(['S Matrix E=',num2str(E),'real psiAIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, imag(psiAIn)); hold on;
title(['S Matrix E=',num2str(E),'imag psiAIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, real(psiBIn)); hold on;
title(['S Matrix E=',num2str(E),'real psiBIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, imag(psiBIn)); hold on;
title(['S Matrix E=',num2str(E),'imag psiBIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, real(psiCIn)); hold on;
title(['S Matrix E=',num2str(E),'real psiCIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, imag(psiCIn)); hold on;
title(['S Matrix E=',num2str(E),'imag psiCIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, probIn); hold on;
title(['S Matrix E=',num2str(E),' probIn'],'interpreter','latex')
shading flat; 
colorbar;
set(gca,'fontsize',25,'FontName','times new roman')

figure()
pcolor(xx, yy, ps); hold on;
title(['S Matrix E=',num2str(E),' sum_3 real^2+img^2'])
shading flat; 
colorbar;
set(gca,'fontsize',15,'FontName','times new roman')