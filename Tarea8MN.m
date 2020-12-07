%Tarea 8

theta=1/4; %Se define el valor de theta, intercambiandolo por 0,1/4,1/2,3/4,1,10,100
frontera=[]; %Se define la frontera de la región A-estable, tomando la expresión (1.6)
for a=0:pi/1000:2*pi % Se define alpha en [0,2pi]
    z=(cos(a)+1i*sin(a)-1)/(1-(theta)*(1-cos(a)-1i*sin(a)));
    frontera=[frontera z];
end

aestable=[]; %Se defina la región A-estable, tomando la expresión (1.5)
for a=0:pi/1000:2*pi
    a_z=(-(1/2))/(1-(1/2)*(theta));
    aestable=[aestable a_z];
end

plot(real(frontera), imag(frontera), 'b', real(aestable), imag(aestable), 'r.');
zoom on; grid on; hold on;
plot( [-5 5],[0 0 ], '-k', 'LineWidth', 1  );
plot( [0 0],[-5 5 ], '-k', 'LineWidth', 1  );
axis square; 
xlabel('Re'); ylabel('Im'); title('Método theta para theta=1/4')
