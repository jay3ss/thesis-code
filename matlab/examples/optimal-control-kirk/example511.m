function example511
% Example 5.1-1 from Kirk's Optimal Control Theory, pg 198
% Modified to have final state free.
    init_x = [0, 0];
    init_p = [0, 0];

    eps = 1e-3;
    options = odeset('RelTol', 1e-4, 'AbsTol',[1e-4 1e-4]);
    t_0 = 0; t_f = 2;
    step = 0.4;
    t_segment = 50;
    t_u = linspace(t_0, t_f, t_segment);   % discretize time
    u = ones(1, t_segment);                % guessed initial control  u=1
    max_iterations = 100;                  % Maximum number of iterations

    % Solve using gradient descent
    for i = 1:max_iterations
        % disp(strcat('iteration: ', num2str(i)));
        % Solve states
        [t_x, X] = ode45(@(t, x) stateEq(t, x, u, t_u), [t_0, t_f], init_x, options);
        x1 = X(:, 1);
        x2 = X(:, 2);

        % Solve costates
        [t_p, P] = ode45(@(t, p) costateEq(t, p), [t_0, t_f], init_p, options);
        p1 = P(:, 1);
        p2 = P(:, 2);

        p1 = interp1(t_p, p1, t_x);
        p2 = interp1(t_p, p2, t_x);

        % size(u)
        dh = dhdu(u, p2, t_x, t_u);
        % size(dh)
        % size(u)
        h_norm = dh'*dh;

        % Find the cost
        J(i, 1) = costFunction(u, t_f, t_u);

        % Should I stay or should I go?
        if h_norm < eps
            J(i, 1)
            break;
        else
            u = adjustControl(dh, t_x, u, t_u, step);
        end
    end

    % Let's do some plotting!
    % plot the state variables & cost for each iteration
    figure(1);
    plot(t_x, X ,'-');
    hold on;
    plot(t_u,u,'r:');legend('x(t)', 'u(t)');
    text(.2,0.08,'x_1(t)');
    text(.25,-.1,'x_2(t)');
    text(.2,.4, 'u(t)');
    s = strcat('Final cost is: J=',num2str(J(end,1)));
    text(.4,1, s);
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
    plot(t_p, P,'-');
    xlabel('Costates');
    ylabel('P');
    legend('p_1(t)', 'p_2(t)')

    if i == max_iterations
        disp('Stopped before required residual is obtained.');
    end

    % Function definitions for system
    function dx = stateEq(t, x, u, t_u)
        u = interp1(t_u, u, t);
        dx = zeros(2, 1);
        dx(1) = x(2);
        dx(2) = -x(2) + u;
    end

    function dp = costateEq(t, p)
        dp = zeros(2, 1);
        dp(1) = 0;
        dp(2) = -p(1) + p(2);
    end

    function cost = costFunction(u, t_f, t_u)
        cost = t_f*0.5*u*u'/length(t_u);
    end

    function dh = dhdu(u, p2, t_x, t_u)
        u = interp1(t_u, u, t_x);
        dh = u + p2;
    end

    function new_u = adjustControl(dh, t_x, u, t_u, step)
        % interploate dH/du
        dh = interp1(t_x, dh, t_u);
        new_u = u - step*dh;
    end
end
