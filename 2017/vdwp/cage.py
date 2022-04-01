import numpy as np

def histogram_1lj(host, eps, sig):
    """
    Make the histogram of energy in every point in a cage.
    """
    # はじめからすべての格子点600x600x600を準備しておき、削る?
    # 遅すぎる。話にならない。
    # やっぱり、CageIntegralのwrapperを作ったほうがはやそうだ。
    # しかし、CageIntegralは高機能すぎて、渡すべき情報が多い。今はもっと簡易にしたいのに。
    # C++で簡易版を作るのもいいかもしれないが、目下の用事には時間がかかりすぎるのであきらめる。
    # 前回、一般化相図上の曲線はどうやって描いたの?
    X = np.linspace(-3.0, 3.0, 600)
    Y = np.linspace(-3.0, 3.0, 600)
    Z = np.linspace(-3.0, 3.0, 600)
    X,Y,Z = np.meshgrid(X,Y,Z)
    X = X.flatten()
    Y = Y.flatten()
    Z = Z.flatten()
    xyz = np.vstack([X,Y,Z]).T
    print(xyz.shape)
    for atom in atoms:
        xyz = xyz[np.linalg.norm(atom-xyz) > sig/2]
    print(xyz.shape)    

mode = None
atoms = None
with open("data/12hedra.nx4a") as f:
    for line in f:
        if line[:5] == "@NX4A":
            mode = "NX4A"
        elif mode == "NX4A":
            if atoms is None:
                natom = int(line)
                atoms = []
            elif natom > 0:
                atoms.append([float(x) for x in line.split()])
                natom -= 1
            else:
                mode = None
                atoms = None
atoms = np.array(atoms)[:,:3]
atoms = atoms[np.linalg.norm(atoms, axis=1) < 6.0]

histogram_1lj(atoms, eps=4.5, sig=3.0)


