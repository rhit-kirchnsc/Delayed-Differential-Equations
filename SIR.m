%Exercise 1 - SIR DDE model

clear
clc
%close all

%parameters
tau=1;                  %time people are sick
tspan=[0 20];               %time
delay = true;           %dde if true, ode if false

%initial conditions
S0=99;                   %initial susceptible
I0=1;                    %initial infected
R0=0;                    %initial recovered
y0=[S0;I0;R0];

if delay
    %history function
    f=@(t) [S0;I0;R0];

    sol=dde23(@SIR_DDE_eqns, tau, f, tspan);
    t_val = sol.x;
    y_val = sol.y;
else
    [a, b] = def_ab;

    [t_val, y_val] = ode45(@(t, y) SIR_ODE_eqns(t, y), tspan, y0);
end

%disp(["",length(y_val)])
%disp(["",length(t_val)])

if (length(t_val) ~= length(y_val))
    y_val_interp = interp1(linspace(t_val(1), t_val(end), size(y_val, 2)), y_val', t_val, 'linear')';
else
    y_val_interp = y_val;
end

%disp(["",length(y_val_interp)])
%disp(["",length(t_val)])

if delay
    figure;
    plot(t_val, y_val_interp(1,:),'-g', 'LineWidth', 2, 'DisplayName', 'Susceptible');
    hold on;
    plot(t_val, y_val_interp(2,:),'-r', 'LineWidth', 2, 'DisplayName', 'Infected');
    plot(t_val, y_val_interp(3,:),'-b', 'LineWidth', 2, 'DisplayName', 'Recovered');
    xlabel('Time');
    title(sprintf('SIR Model, tau = %d', tau));
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
    title(['SIR Model, tau = 0']);
    legend;
    grid on;
    hold off;
end

function dydt = SIR_DDE_eqns(t,y,Z)
    priorI=Z(2,:);
    [a,b]=def_ab;

    dSdt=-a*y(1)*priorI;               %susceptible
    dIdt=a*y(1)*priorI-b*y(2);         %infected
    dRdt=b*y(2);                       %recovered

    dydt=[dSdt;dIdt;dRdt];
end

function dydt = SIR_ODE_eqns(t,y)

    [a,b]=def_ab;

    dSdt=-a*y(1)*y(2);               %susceptible
    dIdt=a*y(1)*y(2)-b*y(2);         %infected
    dRdt=b*y(2);                     %recovered
    dydt=[dSdt;dIdt;dRdt];
end

function [a,b] = def_ab
    a=0.01;                  %infection rate
    b=0.2;                  %recovery rate
end
