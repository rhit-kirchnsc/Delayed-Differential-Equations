clc
clearvars
close all

function y = history(t)
    y = [sin(t),sin(2*t),sin(3*t)]; % from txt bk
end

% Define range of tau values
tau_values = [2,pi,12,13,exp(pi),30]; % example tau values - adjust as needed
tspan = [0 100];

% Pre-allocate solution storage
solutions = cell(length(tau_values), 1);

% Loop through tau values
for i = 1:length(tau_values)
    tau = tau_values(i); % Current tau value
    
    % Define delay functions with current tau
    dely = @(t,y) t - tau;
    delyp = @(t,y) t - tau;
    
    % Define DDE function
    ddefun = @(t,y,ydel,ypdel) DDE_system(t,y,ydel,ypdel);
    
    % Solve for current tau
    solutions{i} = ddensd(ddefun, dely, delyp, @history, tspan);
end

% Evaluate solutions
tn = linspace(0,100,1000);
yn = zeros(3, length(tn), length(tau_values)); % 3 components, time points, tau values

for i = 1:length(tau_values)
    yn(:,:,i) = deval(solutions{i}, tn);
end

% Define the DDE system as a separate function for clarity
function yp = DDE_system(t,y,ydel,ypdel)
    L = 100*[-7 1 2; 3 -9 0; 1 2 -6];
    M = 100*[1 0 -3; -1/2 -1/2 -1; -1/2 -3/2 0];
    N = (1/72)*[-1 5 2; 4 0 3; -2 4 1];
    
    yp = L*y + M*ydel + N*ypdel;
end

%Plotting
figure
for i = 1:length(tau_values)
    subplot(length(tau_values),1,i);
    plot(tn, yn(:,:,i)');
    title(['tau = ' num2str(tau_values(i))]);
    xlabel('t');
    ylabel('y');
    legend('y1','y2','y3');
end