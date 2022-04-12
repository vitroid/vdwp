# Constants
names = ["CS1", "CS2", "TS1", "HS1"]
aliases = {"CS1": "sI", "CS2": "sII", "TS1": "sIII", "HS1": "sIV"}
cage_components = {"CS2": (136., {16: 8., 12: 16.}),  # cage 0 is Nw
                   "CS1": (46., {12: 2., 14: 6.}),
                   "TS1": (172, {12: 10, 14: 16, 15: 4}),
                   "HS1": (80, {12: 6, 14: 4, 15: 4})}
ratios = {"CS1": (1.0, 0.0, 0.0),
          "CS2": (0.0, 1.0, 0.0),
          "HS1": (0.0, 0.0, 1.0),
          "TS1": (23 / 43, 0.0, 20 / 43)}
nma_file = {"CS2": "lattice/C15retry/cat.0.renma.dist",
            "CS1": "lattice/A15retry/cat.0.renma.dist",
            "TS1": "lattice/sigma-udachinx1.006/cat.0.renma.dist",
            "HS1": "lattice/Z-Yang2009/cat.0.renma.dist"}
U_e = {"CS2": -54.7751,
       "CS1": -54.3857,
       "HS1": -54.33,  # MT2011
       "TS1": -54.37,
       "Ih": -55.2238,
       }
