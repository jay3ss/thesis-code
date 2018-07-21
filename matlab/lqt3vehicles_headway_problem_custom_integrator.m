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
X0 = [pos1; 23; pos2; 5; pos3; 15];

% Time
tinitial = 0; % seconds
tfinal = 30; % seconds
dt = 0.1; % seconds
t = tinitial:dt:tfinal;
sim_time = tfinal*1000 + 1;

% K = sym('K%d%d', [6, 6]);
% S = sym('S%d', [6, 1]);
K = zeros(length(Kfinal(:)), length(t));
K(:, length(t)) = Kfinal(:);

S = zeros(length(Sfinal(:)), length(t));
S(:, length(t)) = Sfinal(:);

for i = 1:1:sim_time
    S = integrk4(@sdot, S(), ); 
end

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
%     K = interp1(tk, K, t);
%     K = reshape(K, size(A));
%     
%     S = interp1(ts, S, t);
%     S = S';
    
    U = -(R^-1)*B'*K*X - (R^-1)*B'*S;
    dxdt = A*X + B*U;
    dxdt = dxdt(:);
end