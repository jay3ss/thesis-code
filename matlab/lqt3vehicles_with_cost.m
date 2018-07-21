% Example 5.2-3 from Kirk's Optimal Control Theory
clc
clear all
close all

% Constants
% HDES places the vehicles evenly on the road
VLENGTH = 5; % m (vehicle length)
RLENGTH = 1000; % m (road length)
% VMAX = 31.2928; % m/s (70 mph)
VMAX = 20;
NUM_VEHICLES = 3;
HDES = (RLENGTH - NUM_VEHICLES*VLENGTH)/NUM_VEHICLES;
max_its = 100;%00;
step = 0.001;
epsilon = 1e-4;

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
Q = 2.5e-4*eye(size(A));
% Q(1, 1) = 2.5e-4; Q(3, 3) = 2.5e-4; Q(5, 5) = 2.5e-4;
H = Q;
R = 10*eye(3);

% Problem parameters
r = [HDES; VMAX; HDES; VMAX; HDES; VMAX];

pos3 = RLENGTH;
pos2 = pos3 + 200 - RLENGTH;
pos1 = pos2 + 200;

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
tfinal = 180; % seconds

t_segment = 5e4;
tu = linspace(tinitial, tfinal, t_segment);
u = 1*ones(3, t_segment); % guess the inital control is a step function

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
S = sym('S%d', [6, 1]);

% sim_time = tfinal*1000 + 1;
tic
for j = 1:max_its
%     time_start = tic;
    [tk, K] = ode45(@(t, K)kdot(t, K, A, B, Q, R), [tfinal, tinitial], Kfinal, optionsk);
%     K = fliplr(K);
    [ts, S] = ode45(@(t, S)sdot(t, S, A, B, K, Q, R, r, tk), [tfinal, tinitial], Sfinal, options1);
%     S = fliplr(S);
    [tx, X] = ode45(@(t, X)xdot(t, K, S, X, R, A, B, tk, ts), [tinitial, tfinal], X0, options2);

%     K_new = interp1(tk, K, t_new);
%     K_new = reshape(K_new, size(A));
%     
%     S_new = interp1(ts, S, t_new);
%     S_new = S_new';
% 
%     X_new = interp1(tx, X, t_new);
%     X_new = X_new';

%     u(:,j) = -(R^-1)*B'*K_new*X_new - (R^-1)*B'*S_new;
%     tu(j)=t_new;
    dH = dHdu(tu, u, tx, X, tk, K, ts, S, A, B, R);
    H_norm = normdH(dH);

    J(j, 1) = cost_function(tx, X, tu, u, tfinal, Q, R, r);

    if H_norm < epsilon
        J(j, 1)
        break
    else
        u_old = u;
        u = adjust_control(dH, u, step);
    end

    if mod(j, 10) == 0
       msg = ['Currently on iteration ', num2str(j), ' with cost of ', num2str(J(j, 1))];
       disp(msg);
    end
%     time_elapsed = toc(time_start)
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

figure
plot(ts, S)
xlabel('Time') % x-axis label
ylabel('S(t)') % y-axis label
legend('s_{1}(t)','s_{2}(t)', 's_{3}(t)','s_{4}(t)', 's_{5}(t)','s_{6}(t)')
xlim([0, tfinal+1])
grid on

figure
plot(tx, X(:, 1), tx, X(:, 3), tx, X(:, 5))
xlabel('Time') % x-axis label
ylabel('Headways h_{i}(t)') % y-axis label
legend('h_{1}^*(t)','h_{2}^*(t)','h_{3}^*(t)')
xlim([0, tfinal+1])
grid on

figure
plot(tx, X(:, 2), tx, X(:, 4), tx, X(:, 6))
xlabel('Time') % x-axis label
ylabel('Velocities v_{i}(t)') % y-axis label
legend('v_{1}^*(t)','v_{2}^*(t)','v_{3}^*(t)')
xlim([0, tfinal+1])
grid on

figure
plot(tu, u)
xlabel('Time') % x-axis label
ylabel('u^*(t)') % y-axis label
legend('u_{1}^*(t)', 'u_{2}^*(t)', 'u_{3}^*(t)')
xlim([0, tfinal+1])
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

function dxdt = xdot(t, K, S, X, R, A, B, tk, ts)
    K = interp1(tk, K, t);
    K = reshape(K, size(A));
    
    S = interp1(ts, S, t);
    S = S';

    U = -(R^-1)*B'*K*X - (R^-1)*B'*S;
    dxdt = A*X + B*U;
    dxdt = dxdt(:);
end

function dH = dHdu(tu, u, tx, X, tk, K, ts, S, A, B, R)
    X = interp1(tx, X, tu);
    K = interp1(tk, K, tu);  
    S = interp1(ts, S, tu);
    X = X';
    K = reshape(K, [size(A) length(tu)]);
    S = S';
    len = length(tu);
    for i = 1:len
        dH(:,i) = R*u(:,i) + B'*K(:,:,i)*X(:,i) + B'*S(:,i);
    end
end

function H_norm = normdH(dH)
    dH1 = dH(1, :);
    dH2 = dH(2, :);
    dH3 = dH(3, :);
    H_norm = dH1*dH1' + dH2*dH2' + dH3*dH3';
end

function j = cost_function(tx, X, tu, u, tf, Q, R, r)
    lenx = length(tx);
    X = X';

    track = 0;
    for i = 1:lenx
        track = track + (X(:,i) - r)'*Q*(X(:,i) - r);
    end
%     track = (X' - r)*(X' - r)'*Q;
    lenu = length(tu);
    ctrl = 0;
    for i = 1:lenu
        ctrl = ctrl + u(:,i)'*R*u(:,i);
    end
    
    j = tf/lenx*track + tf/lenu*ctrl;
end

function control = adjust_control(dH, u, step)
    control = u - step*dH;
end