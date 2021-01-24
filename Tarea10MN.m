function Tarea10MN
%Tarea 10
%Se establecen las variables del problema con condiciones iniciales:
x1=-7.5;
x2=7.5;
T=20; %Tiempo final
N=2000; %Malla de dos mil puntos
a=0.2; %Valor para el coeficiente "a"
deltax=(x2-x1)/N;

x=zeros(N+4,1);
for j=-1:N+2 %Se establecen 4 celdas fantasmas, 2 a la izquierda (-1 y 0) y 2 a la derecha (N+1 y N+2) para la condición de frontera
    j_s=j+2;
    x(j_s,1)=x1+deltax*(j-1/2);
end

ui=zeros(N+4,1); %Condición inicial con 4 celdas fantasmas, dos a la izquierda y dos a la derecha para extrapolar la información.
u_exacta=zeros(N+4,1); %Solución exacta

for j=-1:N+2
    j_s=j+2;
    if x(j_s,1) < 0
        ui(j_s,1)=1;
    else
        ui(j_s,1)=0;
    end
end

for j=-1:N+2
    j_s=j+2;
    if x(j_s,1) < a*T
        u_exacta(j_s,1)=1;
    else
        u_exacta(j_s,1)=0;
    end
end

%Se definen las graficas de la condición inicial y la solución exacta:
figure(1)
hold on
plot(x,ui,'r --')
plot(x,u_exacta,'k -')
% stop %Variable no definida para detener el código.

%Se define el método upwind de alta resolución, y se guardan los pasos último y penúltimo para el tiempo final y se actualizan:
%con los minmod como sigue:
function D=minmod(a,b)
if a*b > 0 && a<=b
    D=a;
end
if a+b > 0 && b<=a
    D=b;
end
if a+b <= 0
    D=0;
end
CFL=0.9; %Condición CFL menor a 1, que corresponde a "nu" para ganrantizar estabilidad
dt=CFL*deltax/a; %Valor "delta t" que nos garantiza estabilidad de acuerdo al valor de "a"
nu=a*dt/deltax; %De la definición de "nu" se despeja "dt":
t=0;
u_v=ui;
u_n=ui;

while t < T
    for j=1:N
        j_s=j+2;
        %Se establecen los minmod:
        D_j=minmod(u_v(j_s+1,1)-u_v(j_s,1),u_v(j_s)-u_v(j_s-1,1));
        D_jmas1=minmod(u_v(j_s+2,1)-u_v(j_s+1,1),u_v(j_s+1,1)-u_v(j_s,1));
        D_jmenos1=minmod(u_v(j_s,1)-u_v(j_s-1,1),u_v(j_s-1,1)-u_v(j_s-2,1));
        if a > 0
            u_n(j_s,1)=u_v(j_s,1)-nu*(u_v(j_s,1)-u_v(j_s-1,1))-((1/2)*(nu)*(1-nu)*(D_j-D_jmenos1));
        else
            u_n(j_s,1)=u_v(j_s,1)-nu*(u_v(j_s+1,1)-u_v(j_s,1))+((1/2)*(nu)*(1+nu)*(D_jmas1-D_j));
        end
    end
    
    %Se establecen las condiciones de frontera de Neuman dadas, referidas a
    %las 4 celdas fantasma:
    u_n(1:2,1)=u_n(3,1);
    u_n(N:N+1,1)=u_n(N-1,1);
    u_v=u_n;
    t=t+dt;
end
plot(x,u_n)
end
%Se definen los métodos numéricos de la Tarea 9 anterior:
% one-sided, que viene de izquierda a derecha:
CFL=0.9;
dt=CFL*deltax/a;
nu=a*dt/deltax;
t=0;
u_v=ui;
u_n=ui;
while t < T
    if a < 0
        for j=1:N
            j_s=j+2;
            u_n(j_s,1)=u_v(j_s,1)-nu*(u_v(j_s+1,1)-u_v(j_s,1));
        end
    else
        for j=1:N
            j_s=j+2;
            u_n(j_s,1)=u_v(j_s,1)-nu*(u_v(j_s,1)-u_v(j_s-1,1));
        end
    end
    
    %Se establecen las condiciones de frontera de Neuman dadas, referidas a
    %las 4 celdas fantasma:
    u_n(1:2,1)=u_n(3,1);
    u_n(N+3:N+4,1)=u_n(N+2,1);
    u_v=u_n;
    t=t+dt;
end
plot(x,u_n)

% %Lax-Friedrichs:
CFL=0.9;
dt=CFL*deltax/a;
nu=a*dt/deltax;
t=0;
u_v=ui;
u_n=ui;
while t < T
    for j=1:N
        j_s=j+2;
        u_n(j_s,1)=(1/2)*(u_v(j_s-1,1)+u_v(j_s+1,1))-(nu/2)*(u_v(j_s+1,1)-u_v(j_s-1,1));
    end
    u_v=u_n;
    t=t+dt;
end
plot(x,u_n)

%Lax-Wendroff:

CFL=0.9;
dt=CFL*deltax/a;
nu=a*dt/deltax;
t=0;
u_v=ui;
u_n=ui;
while t < T
    for j=1:N
        j_s=j+2;
        u_n(j_s,1)=u_v(j_s,1)-(nu/2)*(u_v(j_s+1,1)-u_v(j_s-1,1))+(nu^2/2)*(u_v(j_s+1,1)-2*u_v(j_s,1)+u_v(j_s-1,1));
    end
    u_n(1:2,1)=u_n(3,1);
    u_n(N+3:N+4,1)=u_n(N+2,1);
    u_v=u_n;
    t=t+dt;
end
plot(x,u_n)

% Beam-Warming:

CFL=1.9; %en este caso la CFL "nu" está entre 0 y 2, por tanto es menor a 2 en este caso
dt=CFL*deltax/a;
nu=a*dt/deltax;
t=0;
u_v=ui;
u_n=ui;
while t < T
    for j=1:N
        j_s=j+2;
        u_n(j_s,1)=u_v(j_s,1)-(nu/2)*(3*u_v(j_s,1)-4*u_v(j_s-1,1)+u_v(j_s-2,1))+(nu^2/2)*(u_v(j_s,1)-2*u_v(j_s-1,1)+u_v(j_s-2,1));
    end
    u_v=u_n;
    t=t+dt;  
end
plot(x,u_n)

leyenda=legend('Condición inicial','Sol. Exacta','upwind','one-sided','Lax-Friedrichs','Lax-Wendroff','Beam-Warming')
set(leyenda,'Location','West')
end