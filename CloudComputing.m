% Define the DDE system with tasks solved
function dydt = ddefunc(t, y, Z)
    % y(1) = QM, y(2) = Q1, y(3) = Q2, y(4) = Q3
    % y(5) = SM (tasks solved by Master), y(6) = S1, y(7) = S2, y(8) = S3
    % Z(:,1:3) = QM(t-tau_offload1), QM(t-tau_offload2), QM(t-tau_offload3)
    % Z(:,4:6) = Q1(t-tau_return1), Q2(t-tau_return2), Q3(t-tau_return3)
    lambda = 15;       % Task arrival rate
    muM = 3;           % Master processing rate
    mu = [4, 5.5, 4.5];    % Processing rates: mu1, mu2, mu3
    k = [0.3; 0.7; 0.7];  % Offloading rates: k1, k2, k3 (column vector)
    r = [0.2; 0.25; 0.8]; % Return rates: r1, r2, r3 (column vector)
    Qmax = [30; 30; 30];   % Max queue lengths: Q1max, Q2max, Q3max (column vector)

    % Extract current states
    QM = y(1);
    Q = y(2:4);  % Q1, Q2, Q3
    % SM = y(5), S1 = y(6), S2 = y(7), S3 = y(8) are cumulative tasks solved

    % Extract delayed states
    QM_delay = Z(1, 1:3);  % QM(t-tau_offload1), QM(t-tau_offload2), QM(t-tau_offload3)
    Q_delay = Z(2:4, 4:6); % Q1(t-tau_return1), Q2(t-tau_return2), Q3(t-tau_return3)
    Q_delay = diag(Q_delay); % Extract diagonal: Q1(t-tau_return1), etc.

    % Initialize dydt as an 8x1 column vector (4 queues + 4 solved tasks)
    dydt = zeros(8, 1);

    % Compute dQM/dt
    offload_term = sum(k .* QM .* (Q < 0.92 * Qmax));  % Offload only if Qi < Qimax
    return_term = sum(r .* Q_delay .* (Q_delay > 0.90 * Qmax));  % Return if Qi > 0.9*Qimax
    dydt(1) = lambda - muM - offload_term + return_term;

    % Compute dQi/dt for each processor
    for i = 1:3
        offload_in = k(i) * QM_delay(i)  * (Q(i) < Qmax(i));  % Receive from Master 
        return_out = r(i) * Q(i) * (Q(i) > 0.90 * Qmax(i));   % Return to Master
        base_rate = -mu(i) + offload_in - return_out;
        if Q(i)>0
        dydt(i + 1) = (base_rate);  % Prevent decrease below 0
        else
        dydt(i + 1) = max(0, base_rate);
        end
    end
    % Compute rates of tasks solved (dS/dt)
    dydt(5) = muM * all(Q >= 0.99*Qmax);  % Master solves tasks only if QM > 0
    for i = 1:3
        dydt(i + 5) = mu(i) * (Q(i) > 0);  % Pi solves tasks only if Qi > 0
    end
end

% History function (initial conditions for t <= 0)
function y = history(t)
    y = [10; 5; 2; 2; 0; 0; 0; 0];  % QM = 10, Q1 = Q2 = Q3 = 0, SM = S1 = S2 = S3 = 0
end

% Main script
clear variables;
clc;

% Time span
tspan = [0 94];

% Delays: tau_offload (shorter) for QM to Qi, tau_return (longer) for Qi to QM
tau_offload = [1.5, 1.2, 1.3];  % tau_i for offloading
tau_return = [1.7, 2.2, 1];   % tau_j for returns, longer than tau_offload
delays = [tau_offload, tau_return];  % 6 delays total

% Solve the DDE system
sol = ddesd(@ddefunc, delays, @history, tspan);

% Plot queue lengths
figure(1);
subplot(2, 1, 1);
plot(sol.x, sol.y(1,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'QM (Master Queue)');
hold on;
plot(sol.x, sol.y(2,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Q1 (P1 Queue)');
plot(sol.x, sol.y(3,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Q2 (P2 Queue)');
plot(sol.x, sol.y(4,:), 'm-', 'LineWidth', 1.5, 'DisplayName', 'Q3 (P3 Queue)');
hold off;
xlabel('Time (s)');
ylabel('Queue Length');
title('Queue Lengths Over Time');
legend('show');
grid on;

% Plot tasks solved
subplot(2, 1, 2);
plot(sol.x, sol.y(5,:), 'b-', 'LineWidth', 1.5, 'DisplayName', 'SM (Master Tasks Solved)');
hold on;
plot(sol.x, sol.y(6,:), 'r-', 'LineWidth', 1.5, 'DisplayName', 'S1 (P1 Tasks Solved)');
plot(sol.x, sol.y(7,:), 'g-', 'LineWidth', 1.5, 'DisplayName', 'S2 (P2 Tasks Solved)');
plot(sol.x, sol.y(8,:), 'm-', 'LineWidth', 1.5, 'DisplayName', 'S3 (P3 Tasks Solved)');
hold off;
xlabel('Time (s)');
ylabel('Cumulative Tasks Solved');
title('Tasks Solved Over Time');
legend('show');
grid on;