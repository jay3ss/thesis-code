
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

% Symbols for solving system of differential equations
syms x1(t) x2(t) u(t) k11(t) k12(t) k21(t) k22(t) s1(t) s2(t)


% Solve K first
K = [k11, k12, k21, k22];
Kdot = diff(K) == -K*A -A'*K - Q + K*B*(R^-1)*B'*K;
Ksol = dsolve(Kdot);

function dkdt = kdot(t, K, A, B, Q, R)
    K = reshape(K, size(A)); % Convert K from a column vector to a matrix
    dkdt = -K*A - A'*K - Q + K*B*(R^-1)*B'*K;
    dkdt = dkdt(:); % Convert K back to a column vector
end

function dsdt = sdot(t, s, A, B, K, Q, R, r)
    s = reshape(s, size(A)); % this may be a source of future trouble due to size of A
    dsdt = -(A' - K*B*(R^-1)*B')*s + Q*r;
    dsdt = dsdt(:);
end

function dxdt = xdot(t, K, S, X, R, B, T)
    K11 = K(:, 1);
    K12 = K(:, 2);
    K21 = K(:, 3);
    K22 = K(:, 4);
    
    % Interpolation for handling different time steps
    K11 = interp1(T, K11, t);
    K12 = interp1(T, K12, t);
    K21 = interp1(T, K21, t);
    K22 = interp1(T, K22, t);

    K_new = [
        K11, K12;
        K21, K22
    ];

    S1 = S(:, 1);
    S2 = S(:, 2);
    
    S1 = interp1(T, S1, t);
    S2 = interp1(T, S2, t);
    
    S_new = [
        S1;
        S2
    ];

    x1 = X(1);
    x2 = X(2);
    
    % u = -200*(K12*x1 + K22*x2 + s2);
    u = (-R^-1)*B'*K_new - (R^-1)*B'*S_new;
    
    dxdt = zeros(2, 1);
    dxdt(1) = x2;
    dxdt(2) = 2*x1 - x2 + u;
end

