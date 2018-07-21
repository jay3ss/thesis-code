"""Solves the example 5.1-1 from Kirk's Optimal Control Theory
Modified to have final state free.

NOTE: Work in progress!!!
"""
import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt


def x_bc(xa, xb):
    return np.array([xa[0], xb[0]])


def state_eq(t, x, u, tu):
    """Calculates the state equation"""
    # u = u[0, :]
    u_interp = np.interp(t, tu, u[0, :])
    # print(f'u: {u}')
    return np.vstack((x[1], -x[1] + u_interp))
