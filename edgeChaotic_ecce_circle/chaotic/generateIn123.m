function generateIn123(xi,N,xyra)

    R_1 = 1;
    R_2 = 0.6;
    % xi = 0;
    
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
    
    % N=500;
    % xyra = 3;%1.1;
    x_choose=linspace(-xyra,xyra,N);
    y_choose=linspace(-xyra,xyra,N);
    dx = x_choose(2)-x_choose(1);
    
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
    
    save(['In123_xi',num2str(xi),'_N',num2str(N),'_xyra',num2str(xyra),'.mat'],'In1','In2','In3')
end