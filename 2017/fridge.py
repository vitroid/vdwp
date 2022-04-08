# Refigerant data from https://www.chillercity.com/Refrigerant_Data.php

# import requests
# from bs4 import BeautifulSoup
  
# vgm_url = 'https://www.chillercity.com/Refrigerant_Data.php'
# html_text = requests.get(vgm_url).text
# soup = BeautifulSoup(html_text, 'html.parser')
# print(soup)

import pandas as pd
import pickle
from CageIntegral.physconst import NA, NkB


def refrigerants():
    # Cache it
    # url = "https://www.engineeringtoolbox.com/refrigerants-d_902.html"
    # df = pd.read_html(url)
    # with open("fridge.pickle", "wb") as f:
    #     pickle.dump(df, f)

    # Cached data
    with open("fridge.pickle", "rb") as f:
        df = pickle.load(f)
    table = df[0]
    # print(table.columns.values)
    tc = table[('Critical Point', 'Temperature (oF)')]
    pc = table[ ('Critical Point', 'Pressure (psia)')]
    vc = table[('Critical Point', 'Specific Volume (Cu.Ft./lb.)')]
    name0 = table[('Name', 'Name')]
    name = table[('Refrigerant No.', 'Refrigerant No.')]
    mm = table[('Molecular Mass', 'Molecular Mass')]

    # conversion from mad american units
    tc = 5 * (tc-32)/9 + 273.15 # K, sigh
    pc *= 6894.76 # Pa
    vc = vc / 0.453592  # Cu.Ft. / kg
    vc = vc / (0.001 * mm) # Cu.Ft. / mol
    vc = vc * 0.0283168 # m3 / mol
    vc = vc / NA # m3 / molecule
    rhoc = 1 / vc # molecules / m3
    # Scaled Lennard-Jones critical point
    Tc = 1.321 # T*, T*=kT/eps これで温度からepsを算出できる。
    Rhoc = 0.316 # rho*, rho* = rho sig**3
    # 例えば1気圧のメタンの場合。
    # rhoc = 562.2 kg/m3からsigを計算する。 (LJME: sig=3.758e-10 m)
    # print((Rhoc / (562.2 / (16e-3 / NA)))**(1/3))
    # 計算では2.46 AAとなり、だいぶ小さい。うーむ。難しい。

    eps = NkB*tc / Tc # in kJ/mol
    # print(eps)
    Pc = 0.129 # P*, P* = P sig**3 / eps

    sig = (Pc * eps*1e3/NA / pc)**(1/3)  # epsをkJ/molからJ / moleculeに換算。これはそれらしい数値。
    # print(sig*1e10)
    # # MeのTc=-82.5℃
    # tc = 273.15-82.5
    # eps = NkB*tc / Tc
    # print(eps)
    # # Meのpc=45.8 atm
    # pc = 45.8*101326
    # sig = (Pc * eps*1e3/NA / pc)**(1/3)
    # print(sig*1e10)
    # # いつも使っている値 1.2355 3.758 
    # # cpから計算した値 1.1999561070939162 3.8119005871504252
    # # まあまあ良い。これでいこう。

    newtable = pd.concat([name, sig*1e10, eps, tc, pc/100000, name0], axis=1).dropna()
    newtable = pd.concat([sig*1e10, name0], axis=1).dropna()
    return newtable

if __name__ == "__main__":
    refrigerants().to_excel("test.xlsx")
