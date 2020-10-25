% Tarea 3

Nx=100; Ny=100; %Puntos de la malla a graficar en x,y
Lx=1; Ly=1; %Longitud de 1x1 para la celda
Dx=Lx/Nx; Dy=Ly/Ny; %Δx y Δy
T=0.05; %Tiempo final dado
g=9.81; %Constante de la aceleración de la gravedad en la troposfera
H=1; %Altura promedio del agua
r=0.2; %Radio de la perturbación
f=1; %Parámetro de coriolis

%Se definen las condiciones inciales:
x=zeros(Nx+2,Ny+2);
y=zeros(Nx+2,Ny+2);
eta_ini=zeros(Nx+2,Ny+2);
u_ini=zeros(Nx+2,Ny+2);
v_ini=zeros(Nx+2,Ny+2);

for iy=0:Ny+1
    iy_s=iy+1;
    yj=Dy*(iy-1/2);
    for ix=0:Nx+1
        ix_s=ix+1;
        xi=Dx*(ix-1/2);
        x(ix_s,iy_s)=xi;
        y(ix_s,iy_s)=yj;
        eta_ini(ix_s,iy_s)=max(10*(r^2 -(xi-Lx/2)^2 -(yj-Ly/2)^2),0); %El máximo entre "q" y "0"
        %i.e. max(q,0), para este caso q=0 si la parabola concava hacia
        %abajo <= r^2, o q=0 si la norma de dicha parabola > r
    end
end

Dt=0.00001; %Δt cercano a cero para reducir el error

%Se gráfica la condición inicial, una perturbación incial al centro como:
%plot3(x,y,eta_ini)
%stop

%Ahora, se encuentra la solución como sigue:
uE_n=u_ini; uE_nm1=u_ini; %"u" en el punto Este de la celda para "n" y "n+1"
vN_n=v_ini; vN_nm1=v_ini; %"v" en el punto Norte de la celda para "n" y "n+1"
eta_n=eta_ini; eta_nm1=eta_ini; %Los dos pasos para Eta

t=0; %El tiempo que va corriendo
while t<T
    for iy=1:Ny
        iy_s=iy+1;
        
        for ix=1:Nx
            ix_s=ix+1;
        
        %Se define el algoritmo:
        eta_nm1(ix_s,iy_s)=eta_n(ix_s,iy_s)-Dt*H*(uE_n(ix_s,iy_s)-uE_n(ix_s-1,iy_s))/(Dx)-Dt*H*(vN_n(ix_s,iy_s)-vN_n(ix_s,iy_s-1))/(Dy);
        
        uE_nm1(ix_s,iy_s)=uE_n(ix_s,iy_s)+Dt*f*(vN_n(ix_s,iy_s)+vN_n(ix_s,iy_s)+vN_n(ix_s,iy_s-1)+vN_n(ix_s+1,iy_s-1))/(4)-Dt*g*(eta_n(ix_s+1,iy_s)-eta_n(ix_s,iy_s))/(Dx);
        
        vN_nm1(ix_s,iy_s)=vN_n(ix_s,iy_s)-Dt*f*(uE_n(ix_s-1,iy_s)+uE_n(ix_s,iy_s)+uE_n(ix_s-1,iy_s+1)+uE_n(ix_s,iy_s+1))/(4)-Dt*g*(eta_n(ix_s,iy_s+1)-eta_n(ix_s,iy_s))/(Dy);
        end
    end
    
    t=t+Dt;
    eta_n=eta_nm1;
    uE_n=uE_nm1;
    vE_n=vN_nm1;
end
plot3(x,y,H+eta_nm1);