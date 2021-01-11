%Tarea 9
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

%Se definenen los métodos numéricos con:
CFL=0.9; %Condición CFL menor a 1, que corresponde a "nu"
nu=a*dt/deltax; %De la definición de "nu" se despeja "dt":
dt=CFL*deltax/a; %Valor "delta t" que nos garantiza estabilidad de acuerdo al valor de "a"

%Se gráfica la condición inicial y la solución exacta:
figure(1)
hold on
plot(x,ui,'red --')
plot(x,u_exacta,'black -')
t=0;
%Se guardan los pasos último y penúltimo para el tiempo final y se
%actualizan:
u_v=ui;
u_n=ui;

%one-sided, que viene de izquierda a derecha:
% while t < T
%     if a < 0
%         for j=1:N
%             j_s=j+2;
%             u_n(j_s,1)=u_v(j_s,1)-nu*(u_v(j_s+1,1)-u_v(j_s,1));
%         end
%     else
%         for j=1:N
%             j_s=j+2;
%             u_n(j_s,1)=u_v(j_s,1)-nu*(u_v(j_s,1)-u_v(j_s-1,1));
%         end
%     end
%     
%     %Se establecen las condiciones de frontera de Neuman dadas referidas a
%     %las 4 celdas fantasma:
%     u_n(1:2,1)=u_n(3,1);
%     u_n(N+3:N+4,1)=u_n(N+2,1);
%     u_v=u_n;
%     t=t+dt;
% end
% plot(x,u_n)

%Lax-Friedrichs:

% while t < T
%     for j=1:N
%         j_s=j+2;
%         u_n(j_s,1)=(1/2)*(u_v(j_s-1,1)+u_v(j_s+1,1))-nu/2*(u_v(j_s+1,1)-u_v(j_s-1,1));
%     end
%     u_v=u_n;
%     t=t+dt;
% end
% plot(x,u_n)

%Lax-Wendroff:

while t < T
    for j=1:N
        j_s=j+2;
        u_n(j_s,1)=u_v(j_s,1)-nu/2*(u_v(j_s+1,1)-u_v(j_s-1,1))+nu^2/2*(u_v(j+1,1)-2*u_v(j_s,1)+u_v(j_s-1,1));
    end
    u_n(1:2,1)=u_n(3,1);
    u_n(N+3:N+4,1)=u_n(N+2,1);
    u_v=u_n;
    t=t+dt;
end
plot(x,u_n)

%Beam-Warming:
% CFL=0.9; con nu entre 0 y 2, por tanto CFL es menor a 2
% dt=CFL*deltax/a;
% nu=a*dt/deltax;
% t=0;
% u_v=ui;
% u_n=ui;
% while T<t
%     for j=1:N
%         j_s=j+2;
%         u_n(j_s,1)=u_v(j_s,1)-nu/2*(u_v(j_s+1,1)-u_v(j_s-1,1))+nu^2/2*(u_v(j+1,1)-2*u_v(j_s,1)+u_v(j_s-1,1));
%     end
%     u_v=u_n;
%     T=T+dt;
% end
% plot(x,u_n,'Beam-Warming')

leyenda=legend('Condición inicial','Sol. Exacta','one-sided') %,'Lax-Friedrichs','Lax-Wendroff') %,'Beam-Warming')
set(leyenda,'Location','Southwest')
            