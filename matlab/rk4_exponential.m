% Solve y'(t)=-2y(t) with y0=3, 4th order Runge Kutta
y0 = 3;                     % Initial Condition
dt_vals=[0.8, 0.5, 0.1];   % Time step

for j = 1:length(dt_vals)
    dt = dt_vals(j);
    t = 0:dt:2;
    
    ystar = zeros(size(t));  % Preallocate array (good coding practice)
    ystar(1) = y0;           % Initial condition gives solution at t=0.
    
    for i=1:(length(t)-1)
        ystar(i+1) = integrk4(@ydot, ystar(i), t(i), dt); % Approx soln
    end
    yexact = 3*exp(-2*t);     % Exact solution (in general we won't know this)
    plot(t,abs(ystar-yexact), '+-');
    hold on
end

title('Error')
% plot(t,yexact);
grid on

legend('h_1', 'h_2', 'h_3');

function dydt = ydot(y, t)
    dydt = -2*y;
end