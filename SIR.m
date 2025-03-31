%Exercise 1 - SIR DDE model

clear
clc
%close all

%parameters
tau=15;                   %time people are sick
tspan=[0 40];            %time
delay = true;            %dde if true, ode if false

%initial conditions
S0=999;                   %initial susceptible
I0=1;                    %initial infected
R0=0;                   %initial recovered
y0=[S0;I0;R0];

if delay
    %history function
    [a, b] = def_ab;
    f = @(t) [S0 * exp(-a * I0 * t);  % Approximate susceptible history
          I0 * exp(a * S0 * t);   % Approximate infected history
          R0];   
    
    options = ddeset('RelTol',1e-6,'AbsTol',1e-12);
    sol=dde23(@(t,y,Z) SIR_DDE_eqns(t,y,Z), tau, f, tspan, options);
    t_val = sol.x;
    y_val = sol.y;
else

    [t_val, y_val] = ode45(@(t, y) SIR_ODE_eqns(t, y), tspan, y0);
end

if (length(t_val) ~= length(y_val))
    y_val_interp = interp1(linspace(t_val(1), t_val(end), size(y_val, 2)), y_val', t_val, 'linear')';
else
    y_val_interp = y_val;
end


if delay
    figure;
    plot(t_val, y_val_interp(1,:),'-g', 'LineWidth', 2, 'DisplayName', 'Susceptible');
    hold on;
    plot(t_val, y_val_interp(2,:),'-r', 'LineWidth', 2, 'DisplayName', 'Infected');
    plot(t_val, y_val_interp(3,:),'-b', 'LineWidth', 2, 'DisplayName', 'Recovered');
    xlabel('Time');
    title(sprintf('SIR Model, tau = %d, a = %d percent', tau, int16(a*1000)));
    legend;
    grid on;
    hold off;
else
    figure;
    plot(t_val, y_val_interp(:,1),'-g', 'LineWidth', 2, 'DisplayName', 'Susceptible');
    hold on;
    plot(t_val, y_val_interp(:,2),'-r', 'LineWidth', 2, 'DisplayName', 'Infected');
    plot(t_val, y_val_interp(:,3),'-b', 'LineWidth', 2, 'DisplayName', 'Recovered');
    xlabel('Time');
    title(sprintf('SIR Model, tau = 0, a = %d percent', int16(a)));
    legend;
    grid on;
    hold off;
end

function dydt = SIR_DDE_eqns(t,y,Z)

    [a,b]=def_ab;
    c = 0.01;                   % reinfected

    St=y(1);                    %define S(t)
    It=y(2);                    %define I(t)

    Stau=Z(1);                  %define S(tau)
    Itau=Z(2);                  %define I(tau)

    ItHat=a*St*It;              %define infection rate now
    ItauHat=a*Stau*Itau;        %define infection rate at t=tau

    dSdt=-ItHat;          %susceptible
    dIdt=ItHat-ItauHat;         %infected
    dRdt=ItauHat;         %recovered


    dydt=[dSdt;dIdt;dRdt];
end

function dydt = SIR_ODE_eqns(t,y)

    [a,b]=def_ab;
    
    Irate=a*S;

    dSdt=-a*y(1)*y(2);               %susceptible
    dIdt=a*y(1)*y(2)-b*y(2);         %infected
    dRdt=b*y(2);                     %recovered
    dydt=[dSdt;dIdt;dRdt];
end

function [a,b] = def_ab
    a=0.0005;                  %infection rate
    b=0.01;                  %recovery rate
end
