% Example 5.2-3 from Kirk's Optimal Control Theory
clc
clear all
close all
% System parameters
A = [
  0,  1;
  2, -1
];

B = [
  0;
  1
];

% Problem parameters

r = [
  1;
  0
];

Q = [
  2, 0;
  0, 0
];

H = Q;

R = 0.005;

% Boundary conditions
Kfinal = H;
sfinal = -H*r;
X0 = [0; 0];
tinitial = 0; % seconds
tfinal = 15; % seconds

% Symbols for solving system of differential equations
% Declaring K(t) matrix
% k12 and k21 are essentially equal as K(t) is symmetric
syms k11 k12 k21 k22
K=[k11 k12; k21 k22];
% K_dot = -K*A - A'*K - Q + K*B*(R^-1)*B'*K;
[tk, K] = ode45(@(t, K)kdot(t, K, A, B, Q, R), [tfinal, tinitial], Kfinal);

% new_k_sz = [size(A), length(K)];

% K = reshape(K, new_k_sz);

syms s1 s2
S = [s1; s2];
% S_dot = -(A' - K*B*(R^-1)*B')*S + Q*r;
[ts, S] = ode45(@(t, S)sdot(t, S, A, B, K, Q, R, r, tk), [tfinal, tinitial], sfinal);

[tx, X] = ode45(@(t, X)xdot(t, K, S, X, R, A, B, tk, ts), [tinitial, tfinal], X0);

% Plotting time
figure
plot(tk, K)
xlabel('Time') % x-axis label
ylabel('K(t)') % y-axis label
legend('k_{1}(t)','k_{2}(t)','k_{3}(t)','k_{4}(t)')
xlim([0, 16])
grid on

figure
plot(ts, S)
xlabel('Time') % x-axis label
ylabel('S(t)') % y-axis label
legend('s_{1}(t)','s_{2}(t)')
xlim([0, 16])
grid on

for j = 1:1:15001
    K1 = K(:,1);
    K2 = K(:,2);
    K3 = K(:,3);
    K4 = K(:,4);
    t_new=(j-1)/1000;
    K1=interp1(tk,K1,t_new, 'spline');
    K2=interp1(tk,K2,t_new, 'spline');
    K3=interp1(tk,K3,t_new, 'spline');
    K4=interp1(tk,K4,t_new, 'spline');
    K_new=[K1 K2; K3 K4];
    
    s1 = S(:, 1);
    s2 = S(:, 2);
    
    s1 = interp1(ts, s1, t_new, 'spline');
    s2 = interp1(ts, s2, t_new, 'spline');
    
    S_new = [
        s1;
        s2
    ];

    x1 = X(:,1);
    x2 = X(:,2);
    x1=interp1(tx, x1, t_new, 'spline');
    x2=interp1(tx, x2, t_new, 'spline');
    X_new=[x1; x2];

    u(j) = (-R^-1)*(K2*x1 + K4*x2 + s2);
    tu(j)=t_new;
end

figure
plot(tx, X)
xlabel('Time') % x-axis label
ylabel('States u(t)') % y-axis label
legend('x_{1}^*(t)','x_{2}^*(t)')
xlim([0, 16])
grid on

figure
plot(tu, u)
xlabel('Time') % x-axis label
ylabel('u^*(t)') % y-axis label
legend('u^*(t)')
xlim([0, 16])
grid on

function dkdt = kdot(t, K, A, B, Q, R)
    K = reshape(K, size(A)); % Convert K from a column vector to a matrix
    dkdt = -K*A - A'*K - Q + K*B*(R^-1)*B'*K;
    dkdt = dkdt(:); % Convert K back to a column vector
end

function dsdt = sdot(t, s, A, B, K, Q, R, r, tk)
%     s = reshape(s, size(A)); % this may be a source of future trouble due to size of A
    K = interp1(tk, K, t, 'spline');
    K_new = reshape(K, size(A));
%     K1 = K(:, 1);
%     K2 = K(:, 2);
%     K3 = K(:, 3);
%     K4 = K(:, 4);
    
    % Interpolation for handling different time steps
%     K1 = interp1(tk, K1, t);
%     K2 = interp1(tk, K2, t);
%     K3 = interp1(tk, K3, t);
%     K4 = interp1(tk, K4, t);

%     K_new = [
%         K1, K2;
%         K3, K4
%     ];
    dsdt = -(A' - K_new*B*(R^-1)*B')*s + Q*r;
    dsdt = dsdt(:);
end

function dxdt = xdot(t, K, S, X, R, A, B, tk, ts)
    K = interp1(tk, K, t, 'spline');
%     K = K';
    % K1 = K(:, 1);
    K2 = K(:, 2);
    % K3 = K(:, 3);
    K4 = K(:, 4);
    
    % Interpolation for handling different time steps
    % K1 = interp1(tk, K1, t);
    % K2 = interp1(tk, K2, t);
    % K3 = interp1(tk, K3, t);
    % K4 = interp1(tk, K4, t);

    % s1 = S(:, 1);
    s2 = S(:, 2);
    
    % s1 = interp1(ts, s1, t);
    s2 = interp1(ts, s2, t, 'spline');

    x1 = X(1);
    x2 = X(2);
    
    u = (-R^-1)*(K2*x1 + K4*x2 + s2);
    
    dxdt = zeros(2, 1);
    dxdt(1) = x2;
    dxdt(2) = 2*x1 - x2 + u;
end

