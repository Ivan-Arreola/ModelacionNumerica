% Tarea 4

Nx=100; Ny=100; %Puntos de la malla a graficar en x,y
Lx=1; Ly=10; %Longitud de 1x10 para la celda
Dx=Lx/Nx; Dy=Ly/Ny; %Δx y Δy
T=0.25; %Tiempo final
%Valores dados para beta, c y r:
b=1;
c=1;
r=0.2;


%Se definen las condiciones inciales:
u_ini=zeros(Nx+2,Ny+2);
v_ini=zeros(Nx+2,Ny+2);
x=zeros(Nx+2,Ny+2);
y=zeros(Nx+2,Ny+2);
p_ini=zeros(Nx+2,Ny+2);

for iy=0:Ny+1
    iy_s=iy+1;
    yj=-Ly/2+(iy-1/2)*Dy;
    for ix=0:Nx+1
        ix_s=ix+1;
        xi=Dx*(ix-1/2);
        x(ix_s,iy_s)=xi;
        y(ix_s,iy_s)=yj;
    end
end

p_inix=zeros(Nx+2,1);

for ix=0:Nx+1
    ix_s=ix+1;
    xi=x(ix_s,1);
    p_inix(ix_s,1)=max(10*(r^2-(xi-Lx/2)^2),0); %El máximo entre "q" y "0" i.e. max(q,0)
end

for iy=0:Ny+1
    iy_s=iy+1;
    yj=y(1,iy_s);
    for ix=0:Nx+1
        ix_s=ix+1;
        p_ini(ix_s,iy_s)=p_inix(ix_s,1)*exp(-(yj^2)/2);
    end
end

% Se gráfica la condición inicial, una perturbación incial al centro como:
u_ini=p_ini; %SE CAMBIA A u_ini=0 para inciso B)
Dt=1*10^-5; %Δt cercano a cero para reducir el error
% plot3(x,y,p_ini)
% title('Condición inicial')

%Ahora, se encuentra la solución como sigue:
uE_n=u_ini; uE_nm1=uE_n; %"U" en el punto Este de la celda para "n" y "n+1"
vN_n=v_ini; vN_nm1=vN_n; %"V" en el punto Norte de la celda para "n" y "n+1"
p_n=p_ini; p_nm1=p_n; %Los dos pasos para P

t=0; %El tiempo que va corriendo
while t<T
    for iy=1:Ny
        iy_s=iy+1;
        
        for ix=1:Nx
            ix_s=ix+1;
        
        %Se define el algoritmo:
        
        uE_nm1(ix_s,iy_s)=uE_n(ix_s,iy_s)+(b*yj*Dt)*(vN_n(ix_s,iy_s)+vN_n(ix_s,iy_s-1)+vN_n(ix_s+1,iy_s)+vN_n(ix_s+1,iy_s-1))/4-(Dt/Dx)*(p_n(ix_s+1,iy_s)-p_n(ix_s,iy_s));
        
        vN_nm1(ix_s,iy_s)=vN_n(ix_s,iy_s)-((b*yj*Dt)/4)*(uE_n(ix_s-1,iy_s)+uE_n(ix_s,iy_s)+uE_n(ix_s-1,iy_s+1)+uE_n(ix_s,iy_s+1))-(Dt/Dy)*(p_n(ix_s,iy_s+1)-p_n(ix_s,iy_s));
        
        p_nm1(ix_s,iy_s)=p_n(ix_s,iy_s)-(Dt/Dx)*(uE_n(ix_s,iy_s)-uE_n(ix_s-1,iy_s))-(Dt/Dy)*(vN_n(ix_s,iy_s)-vN_n(ix_s,iy_s-1));

        end
    end
   
    p_n=p_nm1;
    uE_n=uE_nm1;
    vN_n=vN_nm1;
    
    t=t+Dt;
end

% u=uE_nm1.*y;
% v=p_nm1.*y;
% p=p_nm1.*y;
% Se gráfica la solución a tiempo final:
plot(x,p_nm1)
% [x,y]=meshgrid(0:0.2:1,0:0.2:10);
% figure
% quiver(x,y,u,v)
% title('Solución a tiempo final')