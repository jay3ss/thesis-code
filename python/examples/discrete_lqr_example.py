"""Example of how to solve a pendulum system using a discrete LQR controller

Taken from

https://github.com/ssloy/tutorials/blob/master/tutorials/pendulum/lqr.py
"""
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt


def discretelqr(A, B, Q, R):
    """
    Solves the discrete time LQR controller
    x[k+1] = A x[k] + B u[k]
    cost = sum(x[k].T*Q*x[k] + u[k].T*R*u[k])
    using the algebraic Riccati equation
    """
    P = np.matrix(scipy.linalg.solve_discrete_are(A, B, Q, R))
    K = np.matrix(scipy.linalg.inv(B.T*P*B + R)*(B.T*P*A))
    return -K


length = 0.21
radius = 0.006
density = 7856 # kg/m^3
m = (2*length) * (radius**2) * (np.pi/4) * density
M = 0.4
dt = 0.02 # s
g = 9.8 # m/s^2

A = np.matrix([
    [1, dt, 0, 0],
    [0, 1, -(3*m*g*dt) / (7*M + 4*m), 0],
    [0, 0, 1, dt],
    [0, 0, (3*g * (m + M)) / (length * (7*M + 4*m)), 1]
    ])

B = np.matrix([
    [0],
    [7*dt / (7*M + 4*m)],
    [0],
    [-3*dt / (length * (7*M + 4*m))]
    ])

print(A, B)

# Weighting matrix for the states in the cost function
Q = np.matrix('1 0 0 0; 0 0.0001 0 0; 0 0 1 0; 0 0 0 .0001')
# Weighting matrix for the control in the cost function
R = np.matrix('0.0005')

K = discretelqr(A, B, Q, R)
print(f'K: {K}')
print('double c[] = {%f, %f, %f, %f};' % (K[0,0], K[0,1], K[0,2], K[0,3]))

nsteps = 250
tspan = np.linspace(0, 2, nsteps, endpoint=True)
xk = np.matrix('0.2; 0; 0.2; 0')

X = []
T = []
U = []

for t in tspan:
    uk = K*xk
    X.append(xk[0, 0])
    T.append(xk[2, 0])
    v = xk[1, 0]
    force = uk[0, 0]
    accel = force / (M + m)
    u = ((1 - .404) * v + dt*accel) / 0.055 / 10
    U.append(u)
    xk = A*xk + B*uk

plt.plot(tspan, X, label='cart position, meters')
plt.plot(tspan, T, label='pendulum angle, radians')
plt.plot(tspan, U, label='control voltage, decavolts')

plt.legend(loc='upper right')
plt.grid()
plt.show()
