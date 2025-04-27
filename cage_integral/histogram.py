# histoをpythonで書いてみる。
import numpy as np


def interaction(sites: np.ndarray, eps: float, sig: float) -> float:
    """
    サイト間相互作用を計算する。

    ゲストは原点にあるものとし、周囲の原子の位置はsitesで与えられる。
    epsとsigはそれぞれゲストと周囲の原子のあいだのLennard-Jones相互作用のパラメータである。
    epsの単位はkJ/mol, sigの単位はÅである。
    周期境界条件はここでは考慮しない。
    """
    d2 = np.sum(sites**2, axis=1)
    d6 = d2**3
    d12 = d6**2
    sig6 = sig**6
    sig12 = sig6**2
    return np.sum(4 * eps * (sig12 / d12 - sig6 / d6))


def scan_inside_a_cage(
    neighbors: np.ndarray, eps: float, sig: float, tick: float
) -> dict:
    """
    ゲストの位置を変えながら、ゲストの相互作用を計算する。

    最初のゲスト位置はケージの中心とし、ゲストと周囲の原子の相互作用を計算し、その値をヒストグラムに追加する。
    相互作用が負の場合は、その位置から上下左右にtickだけずらした位置を順にスキャンする。

    Args:
        sites (np.ndarray): 周囲の原子の位置
        eps (float): ゲストと周囲の原子のあいだのLennard-Jones相互作用のパラメータ
        sig (float): ゲストと周囲の原子のあいだのLennard-Jones相互作用のパラメータ
        tick (float): ゲストの位置をずらす幅

    Returns:
        dict: _description_
    """
    queue = [(0, 0, 0)]
    site_interaction = dict()
    while queue:
        i, j, k = queue.pop()
        if (i, j, k) in site_interaction:
            continue
        offset = np.array([i, j, k]) * tick
        I = interaction(neighbors - offset, eps, sig)
        site_interaction[i, j, k] = I
        if I < 0:
            queue.append((i + 1, j, k))
            queue.append((i, j + 1, k))
            queue.append((i, j, k + 1))
            queue.append((i - 1, j, k))
            queue.append((i, j - 1, k))
            queue.append((i, j, k - 1))
    H = np.histogram(list(site_interaction.values()), bins=1000, range=(-20, 0))
    H = (H[0] * (tick * 1e-10) ** 3, H[1])
    return H


def test():
    import sys
    import json
    import os
    import pathlib

    # 昔のデータ形式は扱いにくいので、このプログラム用にデータを変換する。
    LJME____ = {
        "eps": 1.2355,  # kJ/mol
        "sig": 3.758,  # Å
    }

    TIP4P = {
        "eps": 0.6486943333333333333333333,  # kJ/mol
        "sig": 3.15357794197649532464,  # Å
    }

    # Lorentz-Berthelot rule
    eps = (LJME____["eps"] * TIP4P["eps"]) ** 0.5
    sig = (LJME____["sig"] + TIP4P["sig"]) / 2

    tick = 0.1  # A
    # 0.2Åとしても、J/molの単位までは一致する。

    T = 273.15
    Nk = 0.00831446261815324
    for cage in (12, 14, 16):
        # 水分子の重心位置を酸素の位置とみなしているため、histogram.cppの結果とは微妙に異なる。
        with open(f"{cage}hedra.nx4a", "r") as f:
            neighbors = []
            while True:
                line = f.readline()
                if line == "":
                    break
                columns = line.split()
                if len(columns) == 7:
                    neighbors.append([float(x) for x in columns[:3]])
        neighbors = np.array(neighbors)
        H = scan_inside_a_cage(neighbors, eps, sig, tick)
        I = np.sum(H[0] * np.exp(-(H[1][1:] + H[1][:-1]) / 2 / (Nk * T)))
        f_0 = -Nk * T * np.log(I)
        print(cage, f_0)
        # 単原子分子なら、Pythonで書いてもたいした時間はかからない。
        # しかし、水分子のような多原子分子になると、計算時間が膨大になる。


if __name__ == "__main__":
    test()
