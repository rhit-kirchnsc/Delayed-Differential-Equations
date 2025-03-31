clc
clearvars
close all

global alph beta gamma k

% Read in data
data = readmatrix("NOAA Data.xlsx");
data = flip(data, 1);       % Flip rows
data(1:2, :) = [];            % Remove first row
data_years = data(:, 1);    % Extract years column
data(:, 1) = [];            % Remove years column from data
data = flip(data, 2);       % Flip columns

% Set up constants
 tau1 = [0.01 0.15 0.995 0.9 0.6];
tau2 = [0 0 0 0.1 0.6];
% tau1 = 0.9;
% tau2 = 0.1;
Alpha = [1 1 1 1 1.2];
Beta = [0 0 0 1 0.8];
Gamma = ones(1, 5);
K = [100 100 100 10 10];
taus = [tau1; tau2];

% Define history function with data passed as parameters
history = @(t, y) history1(t, data, data_years);

% History function
function h = history1(t, data, data_years)
    % Assume t in years
    total_months = abs(t * 12);
    additional_months = mod(total_months, 12);
    years = (total_months - additional_months) / 12;
    years = 2023 - years;  % Convert to actual year
    data_idx = years == data_years;  % Find matching year
    month_vec = 1:12;
    additional_months = additional_months+1;
    h = interp1(month_vec, data(data_idx, :)', additional_months);  % Interpolate
end

% Select parameter set (id = 1 for now)
id = 1;
delays = taus(:);
alph = Alpha(id);
beta = Beta(id);
gamma = Gamma(id);
k = K(id);

% Define the delay differential equation
function dydt = ddefun(t,y, Z)
    global alph beta gamma k
    ydelay1 = Z(1);  % First delay term
    ydelay2 = Z(2);  % Second delay term
    dydt = -alph * tanh(k * ydelay1) + beta * tanh(k * ydelay2) + gamma * cos(2 * pi * t);
end

% Solve the DDE
sol = dde23(@ddefun, delays, history, [0.00 15]);

% Optional: Plot the solution
figure;
plot(sol.x, sol.y, 'LineWidth', 2);
xlabel('Time (years)');
ylabel('Sea Surface Temperature (C)');
title('DDE Solution');
grid on;