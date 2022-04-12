import numpy as np
import matplotlib.pyplot as plt


def drawLine(A, B, C, ax=plt, xtick=np.linspace(-0.4, 0.6, 100),
             ytick=np.linspace(-0.4, 0.6, 100)):
    """
    Ax + By + C = 0
    """
    if A == 0:
        X = xtick
        Y = np.zeros_like(X) - C / B
        ax.plot(X, Y)
    elif B == 0:
        Y = ytick
        X = np.zeros_like(Y) - C / A
        ax.plot(X, Y)
    else:
        X = xtick
        Y = (-C - A * X) / B
        ax.plot(X, Y)
