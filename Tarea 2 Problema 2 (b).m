% Problema 2.
% (b).

N=100; L=1;
Dx=L/N; c=1; T=0.2;
x=zeros(N+2);
theta_ini=zeros(N+2);
for ix=0:N+1
    ix_s=ix+1;
    x_i=Dx*(ix-(1/2));
    x(ix_s,1)=x_i;
    
    if 0.25<x_i && x_i<0.75
        theta_ini(ix_s,1)=(x_i-0.25)*(0.75-x_i)
    end
end

theta_n=theta_ini;
theta_nm1=theta_ini;
t=0;
Dt=0.9*(Dx/c); %SE PUEDE CAMBIAR A 0.01 PARA QUE SE VEA BIEN, O REPORTAR EL ANTERIOR, O AMBOS
while t<T
    for ix=1:N
        ix_s=ix+1;
        theta_nm1(ix_s,1)=theta_n(ix_s,1)-Dt*c*(theta_n(ix_s+1,1)-theta_n(ix_s-1,1))/(2*Dx);
    end
    
    t=t+Dt;
    theta_n=theta_nm1;
end
plot(x,theta_ini,x,theta_nm1)
