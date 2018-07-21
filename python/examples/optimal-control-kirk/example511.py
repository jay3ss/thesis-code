"""Solves the example 5.1-1 from Kirk's Optimal Control Theory
Modified to have final state free.

NOTE: Work in progress!!!
"""

import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def example511():
    x_boundary_conditions = [0.05, 0.0]
    init_p = [0.0, 0.0]

    epsilon = 1e-3
    t0 = 0
    tf = 2
    step = 0.4
    t_segments = 50
    tu = np.linspace(t0, tf, t_segments)
    u = np.ones((1, t_segments), dtype='float')
    u_old = np.copy(u)
    dH = 0
    J = []

    # For use of lambda expression during call to solve_ivp
    x_args = [u, tu]

    max_iterations = 100

    for i in range(max_iterations):
        # 1) Solve the states
        # NOTE: taken from https://stackoverflow.com/a/48999583
        x_solution = scipy.integrate.solve_ivp(
            # fun=lambda t, x: state_eq(t, x, u, tu),
            fun=lambda t, x: state_eq(t, x, u, tu),
            t_span=[t0, tf],
            y0=init_x)
        x1, x2 = x_solution.y
        tx = x_solution.t

        # 2) Solve the costates
        p_solution = scipy.integrate.solve_ivp(costate_eq, [tf, t0], init_p)
        p1, p2 = p_solution.y
        tp = p_solution.t

        ji = cost_function(u, tu, tf)
        J.append(ji)

        dH = dHdu(u, p2, tx, tu, tp)
        H_norm = float(dH.dot(dH.T))

        if H_norm < epsilon:
            break
        else:
            print(f'Cost of iteration {i+1} is {ji}')
            u_old = np.copy(u)
            u = adjust_control(dH, tx, u, tu, step)

    plt.figure(figsize=(16, 8), dpi=180)
    plt.plot(tx, x_solution.y[0], label='x1')
    plt.plot(tx, x_solution.y[1], label='x2')
    plt.legend()


    plt.figure(figsize=(16, 8), dpi=180)
    plt.plot(tp, p_solution.y[0], label='p1')
    plt.plot(tp, p_solution.y[1], label='p2')
    plt.legend()

    # print(f'u.shape: {u.shape}')
    # print(f'tu.shape: {tu.shape}')

    plt.figure(figsize=(16, 8), dpi=180)
    plt.plot(tu, u[0], label='u')
    plt.legend()

    # plt.figure(figsize=(16, 8), dpi=180)
    # plt.plot(J)

    plt.show()


def state_eq(t, x, u, tu):
    """Calculates the state equation"""
    # u = u[0, :]
    u_interp = np.interp(t, tu, u[0, :])
    # print(f'u: {u}')
    return [x[1], -x[1] + u_interp]


def costate_eq(t, p):
    """Calculates the costates"""
    return [0.0, -p[0] + p[1]]


def cost_function(u, tu, tf):
    """Calculates the cost"""
    return tf*0.5*u.dot(u.T)/tu.shape[0]


def dHdu(u, p2, tx, tu, tp):
    """Calculates the partial of H (Hamiltonian) w.r.t u"""
    p2 = np.interp(tu, tp, p2)
    return u + p2


def adjust_control(dH, tx, u, tu, step):
    return u - step*dH


if __name__ == '__main__':
    example511()
