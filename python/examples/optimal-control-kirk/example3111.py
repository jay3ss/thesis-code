"""
Example 3.11-1 from Kirk's Optimal Control Theory
Dynamic Programming
"""
import numpy as np
import scipy.linalg as linalg
import matplotlib.pyplot as plt


def riccati(A, B, Q, R):
    """
    Solve the Ricatti equation and return K
    """
    P = np.matrix(linalg.solve_continuous_are(A, B, Q, R))
    return np.matrix(linalg.inv(R)*(B.T*P)) # return K


def control_law(K, X):
    return -2*K*X


# x_dot = x + u
# A = 1
# B = 1
A = np.matrix([[1]])
B = np.matrix([[1]])

Q = np.matrix([[1]])
R = np.matrix([[1]])

K = riccati(A, B, Q, R)
print(f'K: {K[0, 0]}')

nsteps = 2500
tspan = np.linspace(0, 0.0021, nsteps, endpoint=True)
xk = np.matrix([[10]])

X = []
U = []

for t in tspan:
    uk = control_law(K, xk)
    xk = A*xk + B*uk
    U.append(uk[0, 0])
    X.append(xk[0, 0])
    print(f'xk: {xk[0, 0]}  uk: {uk[0, 0]}  t: {t}')

plt.plot(tspan, X, label='State')
plt.plot(tspan, U, label='Control action')

plt.legend()
plt.grid()
plt.show()
