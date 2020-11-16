%Tarea 5

N=100; %Puntos de la malla a graficar
T=1; %Tiempo final dado
k=T/N; %Tamaño de paso dado

%Se definen las condiciones iniciales:
u=zeros(N+1,1); %Aproximación de primer orden (Euler)
u_t=zeros(N+1,1); %Aproximación de segundo orden (trapezoidal)
u_e=zeros(N+1,1); %Solución exacta
t=zeros(N+1,1); %Tiempo que va corriendo

u(1,1)=1;
u_t(1,1)=1;
u_e(1,1)=1;

%Se define el algoritmo con la sol. exacta y los metodos númericos Euler y
%Trapezoidal:
for i=1:N
    ti=(i-1)*k;
    t(i+1,1)=ti+k;
    u(i+1,1)=u(i,1)+k*cos(ti)*u(i,1);
    u_t(i+1,1)=((1+k/2)*cos(ti))/((1-k/2)*cos(ti+k)*u_t(i,1));
    u_e(i+1,1)=exp(sin(ti+k));
end

plot(t,u_e,t,u,t,u_t)
legend('Solución exacta', 'Método de Euler', 'Método Trapezoidal')