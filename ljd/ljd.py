import numpy as np
from vdwp.decorators import deprecated


def cellpotential(
    r: np.ndarray, sigma: float, epsilon: float, z: int, R: float, a: float = 0.0
) -> np.ndarray:
    """
    Lennard-Jones--Devonshire cellpotential function.

    z: coordination number == number of water molecules on the cage surface.
    R: free cavity radius; 3.908, 4.326 AA for 12- and 14-hedral cages, respectively.
    r: position of the guest from the center.
    """

    def delta(N):
        return ((1 - r / R - a / R) ** (-N) - (1 + r / R - a / R) ** (-N)) / N

    ans = (
        2
        * z
        * epsilon
        * (
            sigma**12 / (R**11 * r) * (delta(10) + a * delta(11) / R)
            - sigma**6 / (R**5 * r) * (delta(4) + a * delta(5) / R)
        )
    )
    return ans


def fvalue2(R: float, z: int, sigma: float, epsilon: float, beta: float) -> float:
    """
    Lennard-Jones--Devonshire cellpotential function.

    R: radius
    z: number of water molecules constituting the cavity
    sigma: sigma
    epsilon: epsilon
    beta: beta
    """
    # 0とmaxの間に98点を追加して100点の配列を作る
    r = np.linspace(0, R, 100)
    dr = r[1] - r[0]
    cp = np.zeros(100)
    r1 = r[r < R]
    r1 = r1[1:]
    cp1 = cellpotential(r1, sigma, epsilon, z, R)
    cp[1 : len(cp1) + 1] += cp1
    return (
        -np.log(4 * np.pi * np.trapezoid(np.exp(-beta * cp) * r**2 * 1e-30, dx=dr))
        / beta
    )


# test
@deprecated("fvalue2")
def fvalue(Rz: dict, sigma: float, epsilon: float, beta: float, plt=None) -> float:
    """
    Lennard-Jones--Devonshire cellpotential function.

    Rz: dict of {R: z}; R is the radius and z is the number of water molecules constituting the cavity.
    sigma: sigma
    epsilon: epsilon
    beta: beta
    """
    # 0とmaxの間に98点を追加して100点の配列を作る
    r = np.linspace(0, max(Rz.keys()), 100)
    dr = r[1] - r[0]
    cp = np.zeros(100)
    for R, z in Rz.items():
        r1 = r[r < R]
        length = len(r1)
        r1 = r1[1:]
        cp1 = cellpotential(r1, sigma, epsilon, z, R)
        cp[1 : length + 1] += cp1
    if plt is not None:
        plt.plot(r, cp * beta)
    # plt.plot(r, np.exp(-cellpotential(r, sigma, epsilon, z, R) * beta)*r**2)
    # print(z, np.trapz(np.exp(-beta * cellpotential(r, sigma, epsilon, z, R)) * r**2, dx=r[1]-r[0]))
    return (
        -np.log(4 * np.pi * np.trapezoid(np.exp(-beta * cp) * r**2 * 1e-30, dx=dr))
        / beta
    )


# 2022-03-25とりあえず書いただけ。まだ使っていない。
# 昔、これとは全く別に検証したはずだが、どこにあるかわからない。
# LJDを使ったとしても、数値積分は必要なので、実際効率がいいのかどうかはわからん。
# でも、単原子分子ならリアルタイムでfを計算できるようになる見込み。
# Devonshireの計算は、ケージより外のホスト格子のことは考慮していないのに、けっこうきれいに
# 全原子計算と合うということは、全原子計算もケージまででいいのかもしれないね。
