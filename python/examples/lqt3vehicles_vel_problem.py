"""
Solving an LQT problem with three vehicles
"""
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def sdot(t, S, A, B, K, Q, R, r):
    R_inv = np.linalg.inv(R)
    K = np.reshape(K, A.shape)
    dsdt = -(A.T - K*B*R_inv*B.T)*S + Q*r
    return np.reshape(dsdt, S.shape)


def xdot(t, K, S, X, R, A, B):
    # K_shape = K.shape
    K = np.reshape(K, A.shape)
    R_inv = np.linalg.inv(R)

    U = -R_inv*B.T*K*X - R_inv*B.T*S
    dxdt = A*X + B*U
    return np.reshape(dxdt, X.shape)
