"""Example of how to interpolate a sine function using NumPy"""
import numpy as np
import matplotlib.pyplot as plt


def interpsin():
    # Equivalent call in MATLAB: x = 0:10
    x = np.arange(0, 10, 1)
    y = np.sin(x)

    # Equivalent calls in MATLAB
    # xi = 0:0.25:10
    # yi = interp1(x, y, xi)
    xi = np.arange(0, 10, 0.25)
    yi = np.interp(xi, x, y)

    return [(x, y), (xi, yi)]


def plot(data, interp, func_name):
    x, y = data
    xi, yi = interp
    plt.figure(figsize=(18, 6), dpi=180)
    data, = plt.plot(x, y, 'o', label='Data points')
    interp, = plt.plot(xi, yi, '-x', label='Interpolated values')
    plt.title(f'Interpolation of {func_name}')
    plt.xlabel('x')
    plt.ylabel(f'y = {func_name}')
    plt.legend(handles=[data, interp])
    plt.show()


if __name__ == '__main__':
    data, interp = interpsin()
    plot(data, interp, 'sin(x)')
