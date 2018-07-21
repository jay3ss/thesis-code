clc
clear all
close all

% epsilon = 1e-4;
% epsilon = 1.5;
% epsilon = 2;
epsilon = 1e-3;
% rel_tol = 1e-6;
% abs_tol = [1e-6 1e-6 1e-6 1e-6 1e-6 1e-6];
% ode_options = odeset('RelTol', rel_tol, 'AbsTol', abs_tol);

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

VLENGTH = 5; % m (vehicle length)
R = 70;
RLENGTH = R*2*pi; % m (road length)
VMAX = 5; % m/s (70 mph)
NUM_VEHICLES = 3;
HDES = (RLENGTH - NUM_VEHICLES*VLENGTH)/NUM_VEHICLES;

road_length = 1000; % meters
hsafe = 2 * 2*pi/road_length; % angular safety distance

R = 100;
step = 0.001;
t0 = 0;
tf = 150;
t_segment = 500;
tu = linspace(t0, tf, t_segment);

u = 10*ones(3, t_segment); % guess the inital control is a step function
init_x = [965, 15, 10, 15, 20, 15];
init_p = zeros(1, 6);
max_its = 1000;

for i = 1:max_its
   [tx, X] = ode45(@(t, x)stateEq(t, x, tu, u), [t0, tf], init_x, options2);
   Theta = ring_position(X, road_length);
   [tp, P] = ode45(@(t, p)costateEq(t, p, tx, Theta, hsafe), [tf, t0], init_p, options1);
   P = fliplr(P); % costates are reversed, right? Should they be
   % flipped?
   % x1 = X(:, 1); x2 = X(:, 2);
   % p1 = P(:, 1); p2 = P(:, 2);

   dH = dHdu(tu, u, tp, P(:, 2), R);
   i
   H_norm = normdH(dH)
   % disp(H_norm);

   J(i, 1) = cost(tx, Theta, tu, u, tf, R);

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

% plot the state variables & cost for each iteration
figure(1);
% plot(tx, X ,'-');
% legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'x_4(t)', 'x_5(t)', 'x_6(t)', ...
%        'Location', 'southwest');
plot(tx, X(:, 1), tx, X(:, 3), tx, X(:, 5))
xlabel('Time') % x-axis label
ylabel('Headways h_{i}(t)') % y-axis label
legend('h_{1}^*(t)','h_{2}^*(t)','h_{3}^*(t)')
% xlim([0, length(tx)+1])
grid on
% text(.2,0.08,'x_1(t)');
% text(.25,-.1,'x_2(t)');
% text(.2,.4, 'u(t)');
s = strcat('Final cost is: J=',num2str(J(end,1)));
text(.4,1,s);
xlabel('time');
ylabel('states');
hold off;
print -depsc -r300 eg2_descent.eps

figure(2);
plot(tx, Theta, '-');
hold on;
plot(tu,u, '--');
legend('x_1(t)', 'x_2(t)', 'x_3(t)', 'x_4(t)', 'x_5(t)', 'x_6(t)', ...
       'u_1(t)', 'u_2(t)', 'u_3(t)', 'Location', 'southeast');
text(.4,1,s);
xlabel('time');
ylabel('states (theta)');
hold off;
print -depsc -r300 eg2_descent_theta.eps

figure(3)
plot(J,'x-');
xlabel('Iteration number');
ylabel('J');
print -depsc -r300 eg2_iteration.eps

figure('PaperPositionMode', 'auto');
plot(tp, P,'-');
xlabel('Costates');
ylabel('P');
legend('p_1(t)', 'p_2(t)', 'p_3(t)', 'p_4(t)', 'p_5(t)', 'p_6(t)')
print(gcf, '-dpdf', 'costates.pdf')

% Let's plot the headways and velocities separately from each other
figure(5);
H = zeros(3, length(tx));
H(1, :) = Theta(:, 1);
H(2, :) = Theta(:, 3);
H(3, :) = Theta(:, 5);
plot(tx, H);
xlabel('Time (s)');
ylabel('Headways (angular)');
legend('h_1(t)', 'h_2(t)', 'h_3(t)');

figure(6);
V = zeros(3, length(tx));
V(1, :) = Theta(:, 2);
V(2, :) = Theta(:, 4);
V(3, :) = Theta(:, 6);
plot(tx, V);
xlabel('Time (s)');
ylabel('Velocities');
legend('v_1(t)', 'v_2(t)', 'v_3(t)');


function dx = stateEq(t, x, tu, u)
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
    u = interp1(tu, u', t);
    u = u';
    dx = zeros(6, 1);

    dx(1) = x(6) - x(2);
    dx(2) = u(1);
    dx(3) = x(2) - x(4);
    dx(4) = u(2);
    dx(5) = x(4) - x(6);
    dx(6) = u(3);
    % disp(dx);end
end

function dp =  costateEq(t, p, tx, x, hsafe)
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