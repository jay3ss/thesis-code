"""Solving the pendulum problem using 4th order Runge-Kutta"""
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
import rk4_example as rk4


def derivative(y, t, b=0.25, c=5):
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return np.array(dydt)


def main():
    t0 = 0
    tf = 10

    b = 0.25
    c = 5.0
    yinit = [np.pi - 0.1, 0.0]
    texact = np.linspace(t0, tf, 101)
    yexact = scipy.integrate.odeint(derivative, yinit, texact, args=(b, c))


    plt.figure(figsize=(18,6), dpi=180)
    colors = ['r', 'b', 'g']
    plt.plot(texact, yexact, label='Exact', color='y')
    hvals = [1, 0.5, 0.1]

    smallest_step = min(hvals)
    total_iterations = (tf - t0)/smallest_step

    for color, h in enumerate(hvals):
        t = np.arange(t0, tf+h, h)
        ystar = np.zeros((2, len(t)))

        ystar[0, 0] = np.pi - 0.1
        for i in range(len(t)-1):
            ystar[:,i+1] = rk4.rk4(derivative, ystar[:,i], t[i], h, (b, c))
            percentage = 100*i/total_iterations
            update_msg = f'Iteration {c+1} is {percentage:.2f} % complete.'
            print(update_msg, end='\r')

        plt.plot(
            t, ystar.T, '--o',  label=[f'h={h}'], linewidth=1,
            markersize=3, markerfacecolor='none', color=colors[color])

    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('y(t)')
    plt.title('Numerical solution of y\'(t)=-2y(t) using fourth-order Runge-Kutta')
    plt.grid()
    plt.show()

    def deriv_wrapper(t, y, b, c):
        return derivative(y, t, b, c)

    # Pendulum ODE
    pode = scipy.integrate.ode(deriv_wrapper)
    pode.set_integrator('lsoda',method='bdf', with_jacobian=False)
    pode.set_initial_value(yinit, t0).set_f_params(b, c)

    dt = 0.1
    P = [[3, 0]]
    T = [pode.t]

    while pode.successful() and pode.t < tf:
        T.append(pode.t+dt)
        P.append(pode.integrate(pode.t+dt))
        print(f'{T[-1]} out of {tf} seconds', end='\r')



    plt.figure(figsize=(18, 6), dpi=180)


    for c, h in enumerate(hvals):
        t = np.arange(t0, tf+h, h)
        ystar = np.zeros((2, len(t)))
        ystar[0, 0] = np.pi - 0.1
        for i in range(len(t)-1):
            ystar[:,i+1] = rk4.rk4(derivative, ystar[:,i], t[i], h)
            percentage = 100*i/total_iterations
            update_msg = f'Iteration {c+1} is {percentage:.2f} % complete.'
            print(update_msg, end='\r')

        plt.plot(
            t, ystar.T,'--o',  label=f'h={h}', linewidth=1,
            markersize=3, markerfacecolor='none', color=colors[c])

    plt.plot(T, P, label='ode', color='y')
    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('y(t)')
    plt.title('Numerical solution of y\'(t)=-2y(t) using fourth-order Runge-Kutta')
    plt.grid()
    plt.show()


if __name__ == '__main__':
    main()
