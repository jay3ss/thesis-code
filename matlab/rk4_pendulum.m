clc
clear all
close all
% Solve y'(t)=-2y(t) with y0=3, 4th order Runge Kutta
y0 = 3;                  % Initial Condition
dt=0.05;                 % Time step
t = 0:dt:45;             % t goes from 0 to 10 seconds.

ystar = zeros(2, length(t));  % Preallocate array (good coding practice)


ystar(1) = pi - 0.1;           % Initial condition gives solution at t=0.
ystar(2) = 0;
b = 0.25;
c = 5;
for i=1:(length(t)-1)
    yn = integrk4(@wdot, ystar(:,i), t(i), dt, b, c); % Approx soln
    ystar(1, i+1) = yn(1);
    ystar(2, i+1) = yn(2);
end

plot(t,ystar);
legend('x1', 'x2');
grid on

% Function definitions
function dwdt = wdot(y, t, b, c)
    theta = y(1);
    omega = y(2);
    x1 = omega;
    x2 = -b*omega - c*sin(theta);
    dwdt = [x1; x2];
end