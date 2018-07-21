% Example 5.2-3 from Kirk's Optimal Control Theory
clc
clear all
close all

% Constants
% HDES places the vehicles evenly on the road
VLENGTH = 5; % m (vehicle length)
RLENGTH = 500; % m (road length)
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
Q = 2.5e-5*eye(size(A));
% Q(1, 1) = 1e-4; Q(3, 3) = 1e-4; Q(5, 5) = 1e-4;
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
P0 = zeros(1, 6);
tinitial = 0; % seconds
tfinal = 150; % seconds

t_segment = 500;
tu = linspace(tinitial, tfinal, t_segment);

U = 10*ones(3, t_segment); % guess the inital control is a step function

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
epsilon = 1e-4;
max_its = 1000;

for i = 1:max_its
   [tx, X] = ode45(@(t, X)xdot(t, X, A, B, tu, U), [tinitial, tfinal], X0, options2);
%    Theta = ring_position(X, road_length);
   [tp, P] = ode45(@(t, p)pdot(t, p, tx, X, HDES), [tfinal, tinitial], P0, options1);
%    P = fliplr(P); % costates are reversed, right? Should they be
   % flipped?
   % x1 = X(:, 1); x2 = X(:, 2);
   % p1 = P(:, 1); p2 = P(:, 2);

   dH = dHdu(tu, U, tp, P(:, 2), R);
   H_norm = normdH(dH);
   % disp(H_norm);

%    headways(1,:) = X(1,:);
%    headways(2,:) = X(3,:);
%    headways(3,:) = X(5,:);
   
   J(i, 1) = cost(tx, X, tu, U, tf, R);

   if H_norm < epsilon
       J(i, 1)
       break
   else
       u_old = u;
       u = adjustControl(dH, u_old, step);
   end

   if mod(i, 10) == 0
       msg = ['Currently on iteration ', num2str(i), ' with cost of ', num2str(J(i, 1))];
       disp(msg);
   end
end

if i == max_its
    disp('Stopped before required residual is obtained.');
end


figure('PaperPositionMode', 'auto');
plot(tp, P,'-');
xlabel('Costates');
ylabel('P');
legend('p_1(t)', 'p_2(t)', 'p_3(t)', 'p_4(t)', 'p_5(t)', 'p_6(t)')
print(gcf, '-dpdf', 'costates.pdf')

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

% Let's plot the headways and velocities separately from each other
% figure(5);
% H = zeros(3, length(tx));
% H(1, :) = Theta(:, 1);
% H(2, :) = Theta(:, 3);
% H(3, :) = Theta(:, 5);
% plot(tx, H);
% xlabel('Time (s)');
% ylabel('Headways (angular)');
% legend('h_1(t)', 'h_2(t)', 'h_3(t)');
% 
% figure(6);
% V = zeros(3, length(tx));
% V(1, :) = Theta(:, 2);
% V(2, :) = Theta(:, 4);
% V(3, :) = Theta(:, 6);
% plot(tx, V);
% xlabel('Time (s)');
% ylabel('Velocities');
% legend('v_1(t)', 'v_2(t)', 'v_3(t)');


function dxdt = xdot(t, X, A, B, tu, U)
%     K = interp1(tk, K, t);
%     K = reshape(K, size(A));
% 
% %     U = -(R^-1)*B'*K*X - (R^-1)*B'*S;
%     dxdt = A*X + B*U;
%     dxdt = dxdt(:);
% end
% 
% function dx = stateEq(t, X, A, B, tu, U)
    % Derivative of state vector
    % dx1 = v3 - v1 = x6 - x2
    % dx2 = u1
    % dx3 = v1 - v2 = x2 - x4
    % dx4 = u2
    % dx5 = v2 - v3 = x4 - x3
%     % dx6 = u3
%     u1 = interp1(tu, u(1, :), t);
%     u2 = interp1(tu, u(2, :), t);
%     u3 = interp1(tu, u(3, :), t);
    U = interp1(tu, U', t);
    U = U';
%     dx = zeros(6, 1);
% 
%     dx(1) = x(6) - x(2);
%     dx(2) = u(1);
%     dx(3) = x(2) - x(4);
%     dx(4) = u(2);
%     dx(5) = x(4) - x(6);
%     dx(6) = u(3);
    % disp(dx);end
    dxdt = A*X + B*U;
    dxdt = dxdt(:);
end

function dp =  pdot(t, p, tx, x, hsafe)
    % Derivative of costate vector
    x = interp1(tx, x, t);
    dp = zeros(6, 1);

    dp(1) = x(1) - hsafe;
    dp(2) = -p(1);
    dp(3) = x(3) - hsafe;
    dp(4) = -p(3);
    dp(5) = x(5) - hsafe;
    dp(6) = -p(5);
end


function dH = dHdu(tu, u, tp, p2, R)
    p2 = interp1(tp, p2, tu); 
    dH = R*u;
    dH(2,:) = dH(2,:) + p2;
end

function H_norm = normdH(dH)
    dH1 = dH(1, :);
    dH2 = dH(2, :);
    dH3 = dH(3, :);
    H_norm = dH1*dH1' + dH2*dH2' + dH3*dH3';
end


function j = cost(tx, x, tu, u, tf, R)
    x1 = zeros(3, 1);
    x1(1) = x(6) - x(2);
    x1(2) = x(2) - x(4);
    x1(3) = x(4) - x(6);
    
   
    % disp(strcat('x1^T*x1', num2str(x1'*x1)))
        
    u_cost = u(1, :)*u(1, :)' + u(2, :)*u(2, :)' + u(3, :)*u(3, :)';
    % disp(strcat('u_cost: ', num2str(u_cost)))
    j = tf*(x1'*x1/length(tx) + R*u_cost/length(tu));
    % disp(strcat('j: ', num2str(j)))
end


function control = adjustControl(dH, u, step)
    control = u - step*dH;
end

function x_theta = ring_position(x, road_length)
    [m, n] = size(x);
    x_theta = zeros(m ,n);
    
    x_theta(1, :) = x(1, :) * 2*pi/road_length;
    x_theta(3, :) = x(3, :) * 2*pi/road_length;
    x_theta(5, :) = x(5, :) * 2*pi/road_length;
    for i = 1:n
        if x_theta(1, i) > 2*pi
            x_theta(1, i) = x_theta(1, i) - 2*pi;
        end
        if x_theta(3, i) > 2*pi
            x_theta(3, i) = x_theta(3, i) - 2*pi;
        end
        if x_theta(5, i) > 2*pi
            x_theta(5, i) = x_theta(5, i) - 2*pi;
        end
    end
end