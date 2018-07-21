"""Solves the example 6.2-2 from Kirk's Optimal Control Theory
Modified to have final state free.

NOTE: Work in progress!!!
"""
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def example622():
    eps = 1e-3

    t0 = 0
    tf = 0.78
    tstep = 0.4
    tsegment = 75

    tu = np.linspace(t0, tf, tsegment)
    u = np.ones((1, tsegment), dtype='float')

    initx = [0.05, 0.0]
    initp = [0.0, 0.0]
    R = 0.1

    J = []

    max_iterations = 100

    for i in range(max_iterations):
        print(f'Iteration {i+1}')
        xsol = scipy.integrate.solve_ivp(
            fun=lambda t, x: state_eq(t, x, tu, u),
            t_span=(t0, tf),
            y0=initx
        )
        x1, x2 = xsol.y
        tx = xsol.t
        # Interpolate the states
        x1 = np.interp(tu, tx, x1)
        x2 = np.interp(tu, tx, x2)

        psol = scipy.integrate.solve_ivp(
            fun=lambda t, p: costate_eq(t, p, [x1, x2], tu, u),
            t_span=(tf, t0),
            y0=initp)
        # print('Finished solving costates')
        # print(f'psol.y: {psol.y}')
        p1, p2 = psol.y
        tp = psol.t
        # Interpolate the costates
        p1 = np.interp(tu, tp, p1)
        p2 = np.interp(tu, tp, p2)

        dH = dHdu([x1, x2], [p1, p2], u, R)
        H_norm = float(dH.dot(dH.T))

        # Calculate cost of iteration
        ji = cost_function([x1, x2], u, R, tf, tu)
        J.append(float(ji))

        if H_norm < eps:
            break
        else:
            u_old = u
            u = adjust_control(u, dH, tstep)

    print(f'tu.shape: {tu.shape}')

    plt.figure(figsize=(16, 8), dpi=180)
    plt.plot(tu, x1, label='x1')
    plt.plot(tu, x2, label='x2')
    plt.plot(tu, u[0], label='u1')
    plt.legend()

    plt.figure(figsize=(16, 8), dpi=180)
    plt.plot(tu, p1, label='p1')
    plt.plot(tu, p2, label='p2')
    plt.legend()

    jrange = np.arange(0, len(J))
    plt.figure(figsize=(16, 8), dpi=180)
    plt.plot(jrange, J, label='cost')
    plt.legend()

    plt.show()


def state_eq(t, x, tu, u):
    x1 = x[0]
    x2 = x[1]
    u = u[0, :]
    u = np.interp(t, tu, u)

    dx1 = -2*(x1 + 0.25) + (x2 + 0.5)*np.exp(25*x1/(x1 + 2)) - (x1 + 0.25)*u
    dx2 = 0.5 - x2 - (x2 + 0.5)*np.exp(25*x1/(x1 + 2))

    return [dx1, dx2]


def costate_eq(t, p, x, tu, u):
    p1, p2 = p
    x1, x2 = x
    # print(f'x1: {x1}')
    # print(f'x1.shape: {x1.shape}')
    # print(f'x2: {x2}')
    # print(f'x2.shape: {x2.shape}')
    # print(f'p1: {p1}')
    # print(f'p2: {p2}')
    # x1 = np.interp(t, tx, x1)
    # x2 = np.interp(t, tx, x2)
    # print(f'u: {u[0, :]}')
    u = u[0, :]
    u = np.interp(t, tu, u)
    # The format for dp1 and dp2 matches that of the book (pg 339)
    dp1 = -2*x1 + 2*p1 \
          - p1*(-2*(x1 + 0.25) + (x2 + 0.5)*np.exp(25*x1/(x1 + 2)) \
          - (x1 + 0.25)*u) + p2*(0.5 -x2 \
          - (x2 + 0.5)*np.exp(25*x1/(x1 + 2)))
    # print(f'dp1: {dp1}')
    dp2 = -2*x2 - p1*np.exp(25*x1/(x1 + 2)) \
          + p2*(1 + np.exp(25*x1/(x1 + 2)))
    # print(f'[dp1, dp2]: {[dp1, dp2]}')
    # print(f'dp1.shape: {dp1.shape}')
    return [dp1[0], dp2[0]]


def cost_function(x, u, R, tf, tu):
    x1 , x2 = x
    return tf*(x1.dot(x1.T) + x2.dot(x2.T) + R*u.dot(u.T))/len(tu)


def dHdu(x, p, u, R):
    x1, x2 = x
    p1, p2 = p
    return 2*R*u - p1*(x1 + 0.25)


def adjust_control(u, dH, tstep):
    return u - tstep*dH


if __name__ == '__main__':
    example622()
