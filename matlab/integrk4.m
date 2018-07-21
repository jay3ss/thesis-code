function sol = integrk4(ode, yn, t, dt, varargin)
% Single time step integrator using the fourth order Runge Kutta
    k1 = ode(yn, t+dt, varargin{:});
    y1 = yn + k1*dt/2;

    k2 = ode(y1, t+dt, varargin{:});
    y2 = yn + k2*dt/2;

    k3 = ode(y2, t+dt, varargin{:});
    y3 = yn + k3*dt;

    k4 = ode(y3, t+dt, varargin{:});
    sol = yn + (k1 + 2*k2 + 2*k3 + k4) * dt/6;
end