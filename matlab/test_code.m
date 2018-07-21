function test_code
  % Asgn3_Q7_Descent    Example 2 of optimal control tutorial.
  %    This example is from D.E.Kirk's Optimal control theory: an introduction
  %    example 6.2-2 on page 338-339.
  %    Steepest descent method is used to find the solution.
  eps = 1e-3;
  options = odeset('RelTol', 1e-4, 'AbsTol',[1e-4 1e-4]);
  t0 = 0; tf = 0.78;
  % R = 0.1;
  step = 0.4;
  t_segment = 50;
  Tu = linspace(t0, tf, t_segment);    % discretize time
  u = ones(3,t_segment);               % guessed initial control  u=1
  initx = [5 0 5 0 5 0];                    % initial values for states
  initp = [0 0 0 0 0 0];                       % initial values for costates
  max_iteration = 100;                 % Maximum number of iterations

  hsafe = 2;
  E = zeros(6, 1);                     % holds the error (h - hsafe)
  ri = 0.1;

  % Controller limits
  Umin = -20;
  Umax = 20;
  Ulims = [Umin, Umax];

  num_vehicles = length(initx);
  total_states = 2*num_vehicles;
  road_length = 1000;

  % v = zeros(1, 3);

  for i = 1:max_iteration
    for j = 1:2:total_states
      init_xi = initx(j:j+1);
      init_pi = initp(j:j+1);
      % 1) start with assumed control u and move forward

      % rel_vel = vl-v;
      [Tx,Xj] = ode45(@(t,x) stateEq(t,x,u,Tu), [t0 tf], init_xi, options);

      % 2) Move backward to get the trajectory of costates
      x = Xj(j,1); v = Xj(j,2);
      
      e = x(j) - hsafe;
      [Tp,Pj] = ode45(@(t,p) costateEq(t,p,u,Tu,x,v,Tx,e), [tf t0], init_pi, options);
      p1 = Pj(:,1); p2 = Pj(:,2);
      % P()
      % Important: costate is stored in reverse order. The dimension of
      % costates may also different from dimension of states
      % Use interploate to make sure x and p is aligned along the time axis
      p1 = interp1(Tp,p1,Tx);


      theta(j) = x(j)*2*pi/road_length;
      if theta(j) > 2*pi
         theta(j) = theta(j) - 2*pi;
      end
    end

    s(1) = x(2) - x(1);
    s(2) = x(3) - x(2);
    s(3) = x(1) - x(3);
    if s(1) < 0
        s(1) = s(1) + road_length;
    end

    if s(2) < 0
        s(2) = s(2) + road_length;
    end

    if s(3) < 0
        s(3) = s(3) + road_length;
    end

    X(i,:) = x(:); % position
    V(i,:) = v(:); % velocity
    S(i,:) = s(:); % headway
    % Calculate deltaH with x1(t), x2(t), p1(t), p2(t)
    dH = pH(x1,p1,Tx,u,Tu);
    H_Norm = dH'*dH;

    % NOTE: this will NOT work as is
%     E(i,:) = errorFunction1(s, hsafe);

    % Calculate the cost function
%     e = E(i,:)';
    J(i,1) = costFunction(e, u, tf, Tx, Tu, ri);
        % if dH/du < epslon, exit
    if H_Norm < eps
     % Display final cost
     J(i,1)
     break;
    else
     % adjust control for next iteration
     u_old = u;
     u = AdjControl(dH,Tx,u_old,Tu,step);
   end
  end

  % plot the state variables & cost for each iteration
  figure(1);
  plot(Tx, X ,'-');
  hold on;
  plot(Tu,u,'r:');
  text(.2,0.08,'x_1(t)');
  text(.25,-.1,'x_2(t)');
  text(.2,.4, 'u(t)');
  s = strcat('Final cost is: J=',num2str(J(end,1)));
  text(.4,1,s);
  xlabel('time');
  ylabel('states');
  hold off;
  print -depsc -r300 eg2_descent.eps

  figure(2);
  plot(J,'x-');
  xlabel('Iteration number');
  ylabel('J');
  print -depsc -r300 eg2_iteration.eps

  figure(3);
  plot(Tx, P,'-');
  xlabel('Costates');
  ylabel('P');
  legend('p_1(t)', 'p_2(t)')

  if i == max_iteration
      disp('Stopped before required residual is obtained.');
  end
end

% State equations
function dx = stateEq(t,x,u,Tu)
  dx = zeros(2,1);
  u = interp1(Tu,u,t); % Interploate the control at time t
  dx(1) = x(2);
  dx(2) = u;
end

% Costate equations
function dp = costateEq(t,p,u,Tu,x1,x2,xt, e)
  dp = zeros(2,1);
%   x1 = interp1(xt,x1,t);   % Interploate the state varialbes
%   x2 = interp1(xt,x2,t);
  e = interp1(xt,e,t);
  u = interp1(Tu,u,t);     % Interploate the control
  dp(1) = e;
  dp(2) = -p(1);
end

% Partial derivative of H with respect to u
function dH = pH(x1,p1,tx,u,Tu)
  % interploate the control
  u = interp1(Tu,u,tx);
  R = 0.1;
  dH = 2*R*u - p1.*(x1 + 0.25);
end

% Adjust the control
function u_new = AdjControl(pH,tx,u,tu,step)
  % interploate dH/du
  pH = interp1(tx,pH,tu);
  u_new = u - step*pH;
end

function cost = costFunction(e, u, tf, Tx, Tu, r)
  cost = tf*(e'*e/length(Tx) + r*u'*u/length(Tu));
end

function error = errorFunction1(h, hsafe)
  % This function calculates the error between the headway (h) and the desired
  % headway (hsafe).
  error = h - hsafe;
end

function error = errorFunction2(hl, h)
  error = hl - h;
end

function error = errorFunction3(vl, v)
  error = vl - v;
end
