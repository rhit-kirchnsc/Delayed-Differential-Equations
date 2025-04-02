%Exercise 2 - El-Nino-La-Nina

clear
clc
%close all

%parameters
params.a=1;            
params.b=0;
params.y=1;
params.k=100;
params.tau1=0.01;
params.tau2=1e-9;
tspan=[0 10];

delay = [params.tau1, params.tau2];
%delay = params.tau1;

%history function

f = @(t) 1;   

%options = ddeset('RelTol',1e-6,'AbsTol',1e-12);
sol=dde23(@(t, T, Z)SIR_DDE_eqns(t,T,Z,params), delay, f, tspan);
t_val = sol.x;
y_val = sol.y;

y_val_interp = y_val;


figure;
plot(t_val, y_val_interp(1,:),'-b', 'LineWidth', 2, 'DisplayName', 'Temp');
ylim([-1 1]);
hold on;
xlabel('Time');
title('El-Nino-La-Nina Model Constant History');
legend;
grid on;
hold off;


function dTdt = SIR_DDE_eqns(t,T,Z,params)    

    dTdt= -params.a*tanh(params.k*Z(1))+params.b*tanh(params.k*Z(2))+params.y*cos(2*pi*t);

end
