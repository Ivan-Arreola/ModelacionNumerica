%Tarea10 - VelocidadGeostrofica del problema 2

%Se definen las constantes:
g=9.81; %Valor aproximado de la aceleración de la gravedad
e=2*10^(-4); %Valor de epsilon dado
om=1; %Valor de la rotación

%Se define el dominio:
Nr=100;
Nth=100;
r1=0.8; %Radio interno del anillo
r2=1; %Radio externo del anillo
dr=(r2-r1)/Nr;
r=zeros(Nr,Nth);
theta=zeros(Nr,Nth);
dth=2*pi/Nth;

for ir=1:Nr
    r(ir,:)=r1+dr*(ir-1/2);
end
for ith=1:Nth
    theta(:,ith)=dth*(ith-1/2);
end

x=zeros(Nr,Nth);
y=zeros(Nr,Nth);
for ith=1:Nth
    for ir=1:Nr
        x(ir,ith)=r(ir,ith)*cos(theta(ir,ith));
        y(ir,ith)=r(ir,ith)*sin(theta(ir,ith));
    end
end

%Se define la velocidad geostrofica y su gráfica con el comando quiver:
z=0.1; %Profundidad
ug=zeros(Nr,Nth);
vg=zeros(Nr,Nth);
for ith=1:Nth
    for ir=1:Nr
        xrth=x(ir,ith);
        yrth=y(ir,ith);
        ug(ir,ith)=(e*g)/(2*om)*(yrth/0.2)*(1/r(ir,ith))*(z^2)/2;
        vg(ir,ith)=(e*g)/(2*om)*(yrth/0.2)*(1/r(ir,ith))*(z^2)/2;
    end
end
quiver(x,y,ug,vg,0.1)