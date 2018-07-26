% Example 5.2-3 from Kirk's Optimal Control Theory
clc
clear all
close all

% Constants
% HDES places the vehicles evenly on the road
VLENGTH = 5; % m (vehicle length)
R = 70;
RLENGTH = R*2*pi; % m (road length)
VMAX = 5; % m/s (70 mph)
NUM_VEHICLES = 3;
HDES = (RLENGTH - NUM_VEHICLES*VLENGTH)/NUM_VEHICLES;

% System parameters
A = [
 0, 1, 0,  0, 0, -1;
 0, 0, 0,  0, 0,  0;
 0, 0, 0, -1, 0,  1;
 0, 0, 0,  0, 0,  0;
-1, 0, 1,  0, 0,  0;
 0, 0, 0,  0, 0,  0
];

B = [
0, 0, 0;
1, 0, 0;
0, 0, 0;
0, 1, 0;
0, 0, 0;
0, 0, 1
];

% Problem parameters
Q = 2.5e-5*eye(size(A));
% Q(1, 1) = 1e-4; Q(3, 3) = 1e-4; Q(5, 5) = 1e-4;
Q(2, 2) = 1e-2; Q(4, 4) = 1e-2; Q(4, 4) = 1e-2;
H = Q;
R = 10*eye(3);

% Problem parameters
r = [HDES; VMAX; HDES; VMAX; HDES; VMAX];

pos3 = RLENGTH;
pos2 = pos3 + 125 - RLENGTH;
pos1 = pos2 + 125;

% Boundary conditions
Kfinal = H;
Sfinal = -H*r;

h1_init = pos3 - pos1;
v1_init = 0;
h2_init = abs(pos2 - pos1);
v2_init = 0;
h3_init = abs(pos3 - pos2);
v3_init = 0;

X0 = [h1_init; v1_init; h2_init; v2_init; h3_init; v3_init];
tinitial = 0; % seconds
tfinal = 500; % seconds

% ODE options
rel_tol = 1e-6;
abs_tol = 1e-6*ones(1, 6);
abs_tol_k = 1e-6*ones(1, 36);
non_neg_idx = 1:6;
options1 = odeset('RelTol', rel_tol, 'AbsTol', abs_tol, 'Refine',10);
optionsk = odeset('RelTol', rel_tol, 'AbsTol', abs_tol_k, 'Refine',10);
% This will force the ODE solver to output zero instead of a negative
% number. It takes an array of indices to force to be non-negative
options2 = odeset(options1,'NonNegative',non_neg_idx);

% Symbols for solving system of differential equations
K = sym('K%d%d', [6, 6]);
[tk, K] = ode45(@(t, K)kdot(t, K, A, B, Q, R), [tfinal, tinitial], Kfinal, optionsk);

S = sym('S%d', [6, 1]);
[ts, S] = ode45(@(t, S)sdot(t, S, A, B, K, Q, R, r, tk), [tfinal, tinitial], Sfinal, options1);

[tx, X] = ode45(@(t, X)xdot(t, X, K, S, R, A, B, tk, ts), [tinitial, tfinal], X0, options2);
dt = 1e-1;
sim_time = tfinal*(dt^-1) + 1;
tic
for j = 1:1:sim_time
    t_new=(j-1)*dt;
    K_new = interp1(tk, K, t_new);
    K_new = reshape(K_new, size(A));
    
    S_new = interp1(ts, S, t_new);
    S_new = S_new';

    X_new = interp1(tx, X, t_new);
    X_new = X_new';

    u(:,j) = -(R^-1)*B'*K_new*X_new - (R^-1)*B'*S_new;
   
    % Make sure that we don't go backwards!
    if u(1,j) <= 0 && X_new(2,1) <= 0
        u(1,j) = 0;
    end
    if u(2,j) <= 0 && X_new(4,1) <= 0
        u(2,j) = 0;
    end
    if u(3,j) <= 0 && X_new(6,1) <= 0
        u(3,j) = 0;
    end
%     u(X_new(2,1) <= 0 & u(1,j) <= 0, j) = 0;
%     u(X_new(4,1) <= 0 & u(1,j) <= 0, j) = 0;
%     u(X_new(6,1) <= 0 & u(3,j) <= 0, j) = 0;
    tu(j)=t_new;
end
toc

% Plotting time
% Not plotting Kij for now, because there's 36 elements
% figure
% plot(tk, K)
% xlabel('Time') % x-axis label
% ylabel('K(t)') % y-axis label
% legend('k_{1}(t)','k_{2}(t)','k_{3}(t)','k_{4}(t)')
% xlim([0, tfinal+1])
% grid on

% figure
% plot(ts, S)
% xlabel('Time') % x-axis label
% ylabel('S(t)') % y-axis label
% legend('s_{1}(t)','s_{2}(t)', 's_{3}(t)','s_{4}(t)', 's_{5}(t)','s_{6}(t)')
% xlim([0, tfinal+1])
% grid on
endtime = 250;
figure('DefaultAxesFontSize',16)
p1 = plot(tx, X(:, 1), tx, X(:, 3), tx, X(:, 5));
p1(1).LineWidth = 2;
p1(1).LineStyle = ':';
p1(2).LineWidth = 2;
p1(2).LineStyle = '--';
p1(3).LineWidth = 2;
p1(3).LineStyle = '-.';
title('Headways')
xlabel('Time (s)') % x-axis label
ylabel('Headways h_{i}(t) (m)') % y-axis label
legend('h_{1}^*(t)','h_{2}^*(t)','h_{3}^*(t)')
xlim([0, endtime])
% xlim([0, tfinal+1])
print -depsc -r300 opt_ctrl_headways.eps
grid on

figure('DefaultAxesFontSize',16)
p2 = plot(tx, X(:, 2), tx, X(:, 4), tx, X(:, 6));
p2(1).LineWidth = 2;
p2(1).LineStyle = ':';
p2(2).LineWidth = 2;
p2(2).LineStyle = '--';
p2(3).LineWidth = 2;
p2(3).LineStyle = '-.';
title('Velocities')
xlabel('Time (s)') % x-axis label
ylabel('Velocities v_{i}(t) (m/s)') % y-axis label
legend('v_{1}^*(t)','v_{2}^*(t)','v_{3}^*(t)', 'Location', 'SouthEast')
xlim([0, endtime])
ylim([0, VMAX+1])
% xlim([0, tfinal+1])
print -depsc -r300 opt_ctrl_velocities.eps
grid on

figure('DefaultAxesFontSize',16)
p3 = plot(tu, u(1,:), tu, u(2,:), tu, u(3,:));
p3(1).LineWidth = 2;
p3(1).LineStyle = ':';
p3(2).LineWidth = 2;
p3(2).LineStyle = '--';
p3(3).LineWidth = 2;
p3(3).LineStyle = '-.';
title('Control trajectories')
xlabel('Time (s)') % x-axis label
ylabel('u^*(t) (m/s^2)') % y-axis label
legend('u_{1}^*(t)', 'u_{2}^*(t)', 'u_{3}^*(t)')
xlim([0, endtime])
% xlim([0, tfinal+1])
print -depsc -r300 opt_ctrl_controls.eps
grid on

function dkdt = kdot(t, K, A, B, Q, R)
    K = reshape(K, size(A)); % Convert K from a column vector to a matrix
    dkdt = -K*A - A'*K - Q + K*B*(R^-1)*B'*K;
    dkdt = dkdt(:); % Convert K back to a column vector
end

function dsdt = sdot(t, s, A, B, K, Q, R, r, tk)
    K = interp1(tk, K, t);
    K_new = reshape(K, size(A));
    dsdt = -(A' - K_new*B*(R^-1)*B')*s + Q*r;
    dsdt = dsdt(:);
end

function dxdt = xdot(t, X, K, S, R, A, B, tk, ts)
    K = interp1(tk, K, t);
    K = reshape(K, size(A));
    
    S = interp1(ts, S, t);
    S = S';

    U = -(R^-1)*B'*K*X - (R^-1)*B'*S;
    dxdt = A*X + B*U;
    dxdt = dxdt(:);
end