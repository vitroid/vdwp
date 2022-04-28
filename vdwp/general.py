import numpy as np
import matplotlib.pyplot as plt


def drawLine(A, B, C, style=".", ax=plt, xtick=np.linspace(-0.6, 1.0, 100),
             ytick=np.linspace(-0.6, 1.0, 100)):
    """
    Ax + By + C = 0
    """
    if A == 0:
        X = xtick
        Y = np.zeros_like(X) - C / B
        ax.plot(X, Y, style)
    elif B == 0:
        Y = ytick
        X = np.zeros_like(Y) - C / A
        ax.plot(X, Y, style)
    else:
        X = xtick
        Y = (-C - A * X) / B
        ax.plot(X, Y, style)
