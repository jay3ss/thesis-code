"""
Solves y'(t)=-2y(t) with y0=3 using the 4th order Runge Kutta
Adapted from http://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
"""
import matplotlib.pyplot as plt
import numpy as np


def rk4(func, yn, t, dt, args=None):
    """4th order Runge Kutta for a first-order system

    Args:
        func: function that defines the derivative of the function to be
              integrated
        t:    current time
        yn:   current y-value
        dt:   integration time-step

    Returns:
        yn+1: next y-value
    """
    if args is None:
        k1 = func(yn, t+dt)
        y1 = yn + k1*dt/2

        k2 = func(y1, t+dt)
        y2 = yn + k2*dt/2

        k3 = func(y2, t+dt)
        y3 = yn + k3*dt

        k4 = func(y3, t+dt)
    else:
        k1 = func(yn, t+dt, *args)
        y1 = yn + k1*dt/2

        k2 = func(y1, t+dt, *args)
        y2 = yn + k2*dt/2

        k3 = func(y2, t+dt,*args)
        y3 = yn + k3*dt

        k4 = func(y3, t+dt, *args)
    return yn + (k1 + 2*k2 + 2*k3 + k4) * dt/6


def deriv(y, t):
    """Equation to be integrated
    y'(t)=-2y(t)
    """
    return -2*y


def main():
    t0 = 0
    tf = 2
    texact = np.linspace(t0, tf)
    yexact = 3*np.exp(-2*texact)

    plt.figure(figsize=(18,6), dpi=180)
    plt.plot(texact, yexact, label='Exact')
    hvals = [0.8, 0.2, 0.05]
    smallest_step = min(hvals)
    total_iterations = (tf - t0)/smallest_step
    for c, h in enumerate(hvals):
        t = np.arange(t0, tf+h, h)
        ystar = np.zeros(t.shape)
        ystar[0] = 3
        for i in range(len(t)-1):
            ystar[i+1] = rk4(deriv, ystar[i], t[i], h)
            percentage = 100*i/total_iterations
            update_msg = f'Iteration {c+1} is {percentage:.2f} % complete.'
            print(update_msg, end='\r')

        plt.plot(t, ystar,'--o',  label=f'h={h}', linewidth=1, markersize=3, markerfacecolor='none')

    plt.legend()
    plt.xlabel('time (s)')
    plt.ylabel('y(t)')
    plt.title('Numerical solution of y\'(t)=-2y(t) using fourth-order Runge-Kutta')
    plt.show()


if __name__ == '__main__':
    main()
