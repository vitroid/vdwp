import numpy as np
import matplotlib.pyplot as plt

def cellpotential(r, sigma, epsilon, z, R, a=0.0):
    """
    Lennard-Jones--Devonshire cellpotential function.

    z: coordination number == number of water molecules on the cage surface.
    R: free cavity radius; 3.908, 4.326 AA for 12- and 14-hedral cages, respectively.
    r: position of the guest from the center.
    """
    def delta(N):
        return ((1 - r / R - a / R)**(-N) - (1 + r / R - a / R)**(-N)) / N

    ans = 2 * z * epsilon * (sigma**12 / (R**11 * r) * (delta(10) + a * delta(
        11) / R) - sigma**6 / (R**5 * r) * (delta(4) + a * delta(5) / R))
    return ans

def fvalue(R, sigma, epsilon, beta, z=20):
    r = np.linspace(0,R,100)[1:99]
    plt.plot(r, cellpotential(r, sigma, epsilon, z, R)*beta)
    # plt.plot(r, np.exp(-cellpotential(r, sigma, epsilon, z, R) * beta)*r**2)
    # print(z, np.trapz(np.exp(-beta * cellpotential(r, sigma, epsilon, z, R)) * r**2, dx=r[1]-r[0]))
    return -np.log(4 * np.pi * np.trapz(np.exp(-beta * cellpotential(r, sigma, epsilon, z, R)) * r**2 * 1e-30, dx=r[1]-r[0])) / beta


# 2022-03-25とりあえず書いただけ。まだ使っていない。
# 昔、これとは全く別に検証したはずだが、どこにあるかわからない。
# LJDを使ったとしても、数値積分は必要なので、実際効率がいいのかどうかはわからん。
# でも、単原子分子ならリアルタイムでfを計算できるようになる見込み。
# Devonshireの計算は、ケージより外のホスト格子のことは考慮していないのに、けっこうきれいに
# 全原子計算と合うということは、全原子計算もケージまででいいのかもしれないね。
