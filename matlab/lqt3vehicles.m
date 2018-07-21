% Example 5.2-3 from Kirk's Optimal Control Theory
clc
clear all
close all

% Constants
VLENGTH = 5; % m (vehicle length)
RLENGTH = 1000; % m (road length)
VMAX = 31.2928; % m/s (70 mph)
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
Q = 0.005*eye(size(A));
H = Q;
R = 0.75*eye(3);

% Problem parameters
r = [HDES; VMAX; HDES; VMAX; HDES; VMAX];

% Q = [2, 0;
%      0, 0];

pos3 = 0;
pos2 = pos3 + 25;
pos1 = pos2 + 100;

% Boundary conditions
Kfinal = H;
Sfinal = -H*r;
X0 = [pos1; 23; pos2; 5; pos3; 35];
tinitial = 0; % seconds
tfinal = 30; % seconds

% Symbols for solving system of differential equations
% Declaring K(t) matrix
syms k11 k12 k13 k14 k15 k16 ...
     k21 k22 k23 k24 k25 k26 ...
     k31 k32 k33 k34 k35 k36 ...
     k41 k42 k43 k44 k45 k46 ...
     k51 k52 k53 k54 k55 k56 ...
     k61 k62 k63 k64 k65 k66 ...
 
K = [k11 k12 k13 k14 k15 k16;
     k21 k22 k23 k24 k25 k26;
     k31 k32 k33 k34 k35 k36;
     k41 k42 k43 k44 k45 k46;
     k51 k52 k53 k54 k55 k56;
     k61 k62 k63 k64 k65 k66];

[tk, K] = ode45(@(t, K)kdot(t, K, A, B, Q, R), [tfinal, tinitial], Kfinal);

syms s1 s2 s3 s4 s5 s6
S = [s1; s2; s3; s4; s5; s6];
% S_dot = -(A' - K*B*(R^-1)*B')*S + Q*r;
[ts, S] = ode45(@(t, S)sdot(t, S, A, B, K, Q, R, r, tk), [tfinal, tinitial], Sfinal);

[tx, X] = ode45(@(t, X)xdot(t, K, S, X, R, A, B, tk, ts), [tinitial, tfinal], X0);

% Plotting time
% Not plotting Kij for now, b/c there's 36 of them
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

sim_time = tfinal*1000 + 1;
for j = 1:1:sim_time
    t_new=(j-1)/1000;
    K_new = interp1(tk, K, t_new);
    K_new = reshape(K_new, size(A));
    
    S_new = interp1(ts, S, t_new);
    S_new = S_new';

    X_new = interp1(tx, X, t_new);
    X_new = X_new';

    u(:,j) = -(R^-1)*B'*K_new*X_new - (R^-1)*B'*S_new;
    tu(j)=t_new;
end

figure
plot(tx, X)
xlabel('Time') % x-axis label
ylabel('States u(t)') % y-axis label
legend('h_{1}^*(t)','v_{1}^*(t)', 'h_{2}^*(t)','v_{2}^*(t)', 'h_{3}^*(t)','x_{3}^*(t)')
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

