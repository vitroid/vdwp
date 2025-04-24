from logging import getLogger, DEBUG, basicConfig

# Constants
names = ["CS1", "CS2", "TS1", "HS1"]
aliases = {"CS1": "sI", "CS2": "sII", "TS1": "sIII", "HS1": "sIV"}
cage_components = {
    "CS2": (136.0, {16: 8.0, 12: 16.0}),  # cage 0 is Nw
    "CS1": (46.0, {12: 2.0, 14: 6.0}),
    "TS1": (172, {12: 10, 14: 16, 15: 4}),
    "HS1": (80, {12: 6, 14: 4, 15: 4}),
}
ratios = {
    "CS1": (1.0, 0.0, 0.0),
    "CS2": (0.0, 1.0, 0.0),
    "HS1": (0.0, 0.0, 1.0),
    "TS1": (23 / 43, 0.0, 20 / 43),
}

# 5.0 MPa, 273K, 3rd column (w/o Pauling entropy)
# TIP4P/ICE
# raw data:  FE (w/o Pauling entropy), Nw, volume / (cm^3/mol)

tip4pice_data = {
    "CS1": "   -59.841454832900  1242  22.723050432819",
    "CS2": "   -59.922551474335  1088  23.039907104034",
    "TS1": "   -59.775409212738  1376  22.788576341391",
    "HS1": "   -59.608521883866  1440  22.826573500133",
}
# CS1: 0.7921471658577371 g cm^-3
# CS2: 0.7812531499681448 g cm^-3
# TS1: 0.789869438544369 g cm^-3
# HS1: 0.7885546203373502 g cm^-3

mu_e = dict()
for structure, s in tip4pice_data.items():
    cols = s.split()
    mu_e[structure] = float(cols[0])
    density = 18 / float(cols[2])
    # logger.debug(f"{structure}: {density} g cm^{-3}")


radii = {
    12: 3.988,
    14: 4.331,
    15: 4.527,
    16: 4.587,
}
nmemb = {12: 20, 14: 24, 15: 26, 16: 28}


# compatibility with old system8
# cage_components = {"sII" : (136.,{16:8.,12:16.}), #cage 0 is Nw
#                    "sI"  : (46., {12:2.,14:6.})}
nma_file = {
    "CS2": "lattice/C15retry/cat.0.renma.dist",
    "sI": "lattice/A15retry/cat.0.renma.dist",
}
U_e = {
    "CS2": -54.7751,
    "sI": -54.3857,
    "Ih": -55.2238,
}

# logger.debug(f"cage radii for LJD approximation: {radii}")
