function simpleode45

[tv, Yv] = ode45(@funsys, [0, pi/2], [1;-1;0]);
% plot(tv, Yv(:,1), '+', tv, Yv(:,2), 'x', tv, Yv(:,3), 'o');
plot(tv, Yv);
hold
grid
title('Simple ODE45')
% text(0.3, 14, '-+- y_1')
% text(0.3, 10, '-+- y_2')
% text(0.3, -12, '-+- y_3')
xlabel('time')
legend('y_1', 'y_2', 'y_3', 'Location', 'northwest')
hold off

function Fv = funsys(t, Y)
  Fv(1,1) = 2*Y(1) + Y(2) + 5*Y(3) + exp(-2*t);
  Fv(2,1) = -3*Y(1) - 2*Y(2) - 8*Y(3) + 2*exp(-2*t) - cos(3*t);
  Fv(3,1) = 3*Y(1) + 3*Y(2) + 2*Y(3) + cos(3*t);
