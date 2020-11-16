%Tarea 5

N=100000; %Puntos de la malla a graficar
T=100; %Tiempo final dado
k=T/N; %Tamaño de paso dado

%Se definen las condiciones iniciales:
u=zeros(N+1,1); %Aproximación de primer orden (Euler)
u_t=zeros(N+1,1); %Aproximación de segundo orden (trapezoidal)
u_e=zeros(N+1,1); %Solución exacta
t=zeros(N+1,1); %Tiempo que va corriendo

u(1,1)=1;
u_t(1,1)=1;
u_e(1,1)=1;

%Se define el algoritmo con la solución exacta y los metodos númericos Euler y
%Trapezoidal:
for n=1:N
    tn=(n-1)*k;
    t(n+1,1)=tn+k;
    u(n+1,1)=u(n,1)+k*cos(tn)*u(n,1);
    u_t(n+1,1)=[(1+k/2*cos(tn))/(1-k/2*cos(tn+k))]*u_t(n,1);
    u_e(n+1,1)=exp(sin(tn+k));
end

plot(t,u_e,'r -',t,u,'b .',t,u_t,'black .')
legend('Solución Exacta','Método de Euler','Método Trapezoidal')
title('Comparación de los métodos de Euler y trapezoidal con la solución exacta a T=100 y N=100000')