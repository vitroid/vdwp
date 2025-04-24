import numpy as np
import matplotlib.pyplot as plt
from logging import getLogger


def drawLine(
    A,
    B,
    C,
    style=".",
    ax=plt,
    xtick=np.linspace(-0.6, 1.0, 100),
    ytick=np.linspace(-0.6, 1.0, 100),
):
    """
    Ax + By + C = 0
    """
    if A == 0 and B == 0:
        # no coexistence
        return
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


def derivatives(x, y):
    logger = getLogger(__name__)
    if len(x) != len(y):
        logger.warning(f"List size is different.")
    n = len(y)
    # intermediate values
    d = np.zeros(n)
    d[1 : n - 1] = (y[2:n] - y[0 : n - 2]) / (x[2:n] - x[0 : n - 2])
    d[0] = (y[1] - y[0]) / (x[1] - x[0])
    d[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    return d


# assume the function increases monotonically
# これはnp.interpolateで置き換えられる。
def zerocross(X, Y):
    if Y[0] > 0:
        return X[0] - 1
    if Y[-1] < 0:
        return X[-1] + 1
    for i in range(0, len(Y) - 1):
        z = Y[i] * Y[i + 1]
        if z <= 0:
            delta = X[i + 1] - X[i]
            slope = (Y[i + 1] - Y[i]) / delta
            dx = -Y[i] / slope
            return X[i] + dx
    print("ERROR zerocross", X, Y)


# find the value at x
def value(x, X, Y):
    return zerocross(Y, X - x)


# X = numpy.array([1,2,3,4,5])
# Y = numpy.array([-3,-1,-4,-1,-5])
# print(value(2.5,X,Y))
