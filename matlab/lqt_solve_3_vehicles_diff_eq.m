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

% Trajectory to track
% r = [h1; v1; h2; v2; h3; v3];
r = [HDES; VMAX; 2*HDES; VMAX; 3*HDES; VMAX]; % evenly spread out vehicles and have them travel at max speed

% syms h1(t) v1(t) h2(t) v2(t) h3(t) v3(t) u1(t) u2(t) u3(t)
syms h1 v1 h2 v2 h3 v3 u1 u2 u3 ...
     s1 s2 s3 s4 s5 s6 ...
     k11 k12 k13 k14 k21 k22 k23 k24 k31 k32 k33 k34 k41 k42 k43 k44

K = [
k11, k12, k13, k14;
k21, k22, k23, k24;
k31, k32, k33, k34;
k41, k42, k43, k44,
];

S = [s1; s2; s3; s4; s5; s6];

% X = [h1; v1; h2; v2; h3; v3];
% U = [u1; u2; u3];
% Xdot = A*X + B*U;

% h1dot = diff(h1) == v1 - v3;
% v1dot = diff(v1) == u1;
% h2dot = diff(h2) == v3 - v2;
% v2dot = diff(v2) == u2;
% h3dot = diff(h3) == v2 - v1;
% v3dot = diff(v3) == u3;
% 
% Xdot = [h1dot; v1dot; h2dot; v2dot; h3dot; v3dot];
% sol = dsolve(Xdot);
% 
% v1(t) = sol.v1;

function dkdt = kdot(t, K, A, B, Q, R)
    K = reshape(K, size(A)); % Convert K from a column vector to a matrix
    dkdt = -K*A - A'*K - Q + K*B*(R^-1)*B'*K;
    dkdt = dkdt(:); % Convert K back to a column vector
end

function dsdt = sdot(t, s, A, B, K, Q, R, r, tk)
    % s = reshape(s, size(A)); % this may be a source of future trouble due to size of A
    K = interp1(tk, K', t);
%     K = K';
    K = reshape(K, size(A));
%     K1 = K(:, 1);
%     K2 = K(:, 2);
%     K3 = K(:, 3);
%     K4 = K(:, 4);
%     K5 = K(:, 5);
%     K6 = K(:, 6);
%     K7 = K(:, 7);
%     K8 = K(:, 8);
%     K9 = K(:, 9);
%     K10 = K(:, 10);
%     K11 = K(:, 11);
%     K12 = K(:, 12);
%     K13 = K(:, 13);
%     K14 = K(:, 14);
%     K15 = K(:, 15);
%     K16 = K(:, 16);
    
    % Interpolation for handling different time steps
%     K1 = interp1(tk, K1, t);
%     K2 = interp1(tk, K2, t);
%     K3 = interp1(tk, K3, t);
%     K4 = interp1(tk, K4, t);
%     K5 = interp1(tk, K5, t);
%     K6 = interp1(tk, K6, t);
%     K7 = interp1(tk, K7, t);
%     K8 = interp1(tk, K8, t);
%     K9 = interp1(tk, K9, t);
%     K10 = interp1(tk, K10, t);
%     K11 = interp1(tk, K11, t);
%     K12 = interp1(tk, K12, t);
%     K13 = interp1(tk, K13, t);
%     K14 = interp1(tk, K14, t);
%     K15 = interp1(tk, K15, t);
%     K16 = interp1(tk, K16, t);

%     K_new = [
%         K1, K2, K3, K4;
%         K5, K6, K7, K8;
%         K9, K10, K11, K12;
%         K13, K14, K15, K16;
%     ];
    dsdt = -(A' - K_new*B*(R^-1)*B')*s + Q*r;
    dsdt = dsdt(:);
end

function dxdt = xdot(t, K, S, X, R, B, tk, ts)
    K = interp1(tk, K, t);
%     K = K';
    K = reshape(K, size(A));
%     K1 = K(:, 1);
%     K2 = K(:, 2);
%     K3 = K(:, 3);
%     K4 = K(:, 4);
%     K5 = K(:, 5);
%     K6 = K(:, 6);
%     K7 = K(:, 7);
%     K8 = K(:, 8);
%     K9 = K(:, 9);
%     K10 = K(:, 10);
%     K11 = K(:, 11);
%     K12 = K(:, 12);
%     K13 = K(:, 13);
%     K14 = K(:, 14);
%     K15 = K(:, 15);
%     K16 = K(:, 16);
    
    % Interpolation for handling different time steps
%     K1 = interp1(tk, K1, t);
%     K2 = interp1(tk, K2, t);
%     K3 = interp1(tk, K3, t);
%     K4 = interp1(tk, K4, t);
%     K5 = interp1(tk, K5, t);
%     K6 = interp1(tk, K6, t);
%     K7 = interp1(tk, K7, t);
%     K8 = interp1(tk, K8, t);
%     K9 = interp1(tk, K9, t);
%     K10 = interp1(tk, K10, t);
%     K11 = interp1(tk, K11, t);
%     K12 = interp1(tk, K12, t);
%     K13 = interp1(tk, K13, t);
%     K14 = interp1(tk, K14, t);
%     K15 = interp1(tk, K15, t);
%     K16 = interp1(tk, K16, t);

%     K_new = [
%         K1, K2, K3, K4;
%         K5, K6, K7, K8;
%         K9, K10, K11, K12;
%         K13, K14, K15, K16;
%     ];

    
end