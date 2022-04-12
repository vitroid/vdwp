#!/usr/bin/env python
#coding: utf-8

import histo2f
import histo
import derivative
import physconst as pc

#とりあえず、pretty printingはあとまわしにして、すべての計算をpython上で行えるようにする。
temperature = 273.15 # Kelvin
pressure     = 101325 # Pa

#! Test case for methane hydrate
import molecule
file = open("data/DEFR")
moldict = dict()
while True:
    line = file.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if len(columns) > 0 and columns[0] in ("@DEFR", "@DEFP"):
        tag = columns[0]
        mol = molecule.Molecule()
        id = mol.load(file, tag)
        moldict[id] = mol


guest = {
    "histo14file":"data/LJME____.14hedra.histo",
    "mass": 16.0,  #molar mass in g/mol
    "symmetry":  1,
    "dimension": 0,        #0 for monatomic, 1 for rodlike, 3 for rigid rotor
    "degreeOfFreedom": 3,  #3 for monatomic, 5 for rodlike, 6 for rigid rotor
}

lattice = {
    "cages":[2,6,0,0],
    "Nw":46,
}



#
#f, mu: free energy/chemical potential
#h    : enthalpy
#
#_w   : of water
#_i   : of ice
#_e   : of the empty lattice 
#_g   : of the guest gas
#_c   : in the cage
#



import numpy

temperatures = numpy.array(range(263,283,2))

#test for enthalpy
#1. fを温度で微分してSを求め、TSを加算することでエンタルピーを求める。

histo14 = histo.loadAHisto( open(guest["histo14file"]) )
#free energy of cage occupation
f14s = numpy.array([histo2f.fvalue(histo14,T) for T in temperatures])
#derivative by T == entropy
S14s = numpy.array(derivative.derivatives(temperatures,f14s))
#Energy of cage occupation
u14s = f14s - temperatures*S14s

#2. ケージ内平均ポテンシャルを計算し、それに圧力項を加算する。
w14s = numpy.array([histo2f.energy(histo14,T) for T in temperatures])
#ここで加算する項は分子の形に依存する。
h2 = w14s + (guest["degreeOfFreedom"]/2.0)*pc.NkB*temperatures

#u14sとw14sは一致するが、h2はそれより大きくなる。正しいのはたぶんh2のほう。
#fの計算はもともと定積の分配関数から導出しているので、それを温度で微分したものには圧力項が含まれない。

import chempot
#mu_g^{id}: chemical potential of the ideal gas
mugids   = numpy.array([chempot.Monatomic( T, pressure, guest["mass"] ) for T in temperatures])

Sg = numpy.array(derivative.derivatives(temperatures,mugids))
hg = mugids - temperatures * Sg

#chemical potential of water ice
#Q構造のポテンシャルと振動の寄与から算出する。
import vdWP
mui = [ -55.2238 + vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", T) for T in temperatures ]
hi = mui - temperatures*derivative.derivatives(temperatures,mui)

#chemical potential of liquid water
#Antoineの式を使う。
import antoine
muw = [ antoine.ChemicalPotentialOfWater(T) for T in temperatures ]
#運動エネルギーの項を含む場合。上とどうちがうのか。
A = 5.40221
B = 1838.675
C = -31.737
tip4p = moldict["TIP4P   "]
muw2 = [chempot.RigidRotor(T,antoine.VaporPressure(T, A,B,C), tip4p.moi[0], tip4p.moi[1], tip4p.moi[2], tip4p.mass, 2.0) for T in temperatures]
hw = muw2 - temperatures*derivative.derivatives(temperatures,muw2)

#水の蒸気圧を確認
print u"水の蒸気圧"
print "T/C, T/K, VP/Pa"
for T in range(0,110,10):
    K = T + 273.15
    print T, K, antoine.VaporPressure(K, A,B,C)

#delta mu
#cPentaine
A = 4.24714
B = 1235.305
C = -30.666
print u"シクロペンタンの蒸気圧"
print "T/C, T/K, VP/Pa"
for T in range(0,110,10):
    K = T + 273.15
    print T, K, antoine.VaporPressure(K, A,B,C)
#value for structure II
muc0_g = [ -54.7751 + vdWP.FreeEnergyOfVibration("data/C15retry/cat.0.renma.dist", T) for T in temperatures ]
guest = moldict["CPENTANE"]
histo16 = histo.loadAHisto( open("data/CPENTANE.16hedra.histo") )
#lUse external pressure
mu_g = [chempot.RigidRotor(T, antoine.VaporPressure(T,A,B,C), guest.moi[0], guest.moi[1], guest.moi[2], guest.mass, 2.0) for T in temperatures]
#mu_g = [chempot.RigidRotor(T,antoine.VaporPressure(T, A,B,C), guest.moi[0], guest.moi[1], guest.moi[2], guest.mass, 2.0) for T in temperatures]
f    = [histo2f.fvalue(histo16,T)+chempot.Avalue(T,guest.mass)+chempot.Cvalue(T,guest.moi[0], guest.moi[1], guest.moi[2])+chempot.SymmetryFix(T,2.0)+chempot.IntegrationFix(T,3) for T in temperatures]

import math
Deltamu=[]
for i in range(len(temperatures)):
    T = temperatures[i]
    beta = 1.0/(pc.NkB*T)
    Deltamu.append(-pc.NkB*T*(8.0/136.0)*math.log(1+math.exp(beta*(mu_g[i]-f[i]))))
Deltamu = numpy.array(Deltamu)
print mu_g
print f
print muc0_g
print Deltamu
mu_g = muc0_g + Deltamu
print mu_g

#結局、実験の蒸気圧をつかった数字がのきなみ信用できない。
#それはさておき、潜熱の計算にいよいよとりかかる。

#空のケージの潜熱
hc0_g = muc0_g - temperatures* numpy.array(derivative.derivatives(temperatures,muc0_g))
#占有率。計算するまでもなく1
y = []
for i in range(len(temperatures)):
    T = temperatures[i]
    beta = 1.0/(pc.NkB*T)
    e = math.exp(beta*(mu_g[i]-f[i]))
    y.append(e/(1.0+e))
y = numpy.array(y)
print y

h_g = mu_g - temperatures* numpy.array(derivative.derivatives(temperatures,mu_g))
w_g = f - temperatures* numpy.array(derivative.derivatives(temperatures,f)) +  (6.0/2.0)*pc.NkB*temperatures

DeltaH = []
for i in range(len(temperatures)):
    T = temperatures[i]
    beta = 1.0/(pc.NkB*T)
    DeltaH.append((8./136.)*y[i]*(h_g[i] - w_g[i]))
DeltaH = numpy.array(DeltaH)
h = hc0_g - DeltaH
    






def test(id08):
    mol = moldict[id08]
    T = 273.15
    print chempot.RigidRotor(T,pressure, mol.moi[0], mol.moi[1], mol.moi[2], mol.mass, 2.0) 
    print id08
    print chempot.Avalue(T,mol.mass) , "A"
    print chempot.Cvalue(T,mol.moi[0], mol.moi[1], mol.moi[2]) , "C"
    print chempot.SymmetryFix(T,2.0), "symm fix"
    print chempot.IntegrationFix(T,3), "Integ fix"

test("CJTHF___")
test("CPENTANE")


import matplotlib.pyplot as plt
plt.plot(temperatures,f14s, label="f_14")
plt.plot(temperatures,u14s, label="u_14")
plt.plot(temperatures,w14s, label="w_14")
plt.plot(temperatures,h2, label="h2_14")
plt.plot(temperatures,hg, label="h_g")
plt.plot(temperatures,mugids, label="mu_g^{id}")
plt.title("Test(Methane)")
plt.xlabel("Temperature [K]")
plt.ylabel("F | U | H  [kJ/mol]")
plt.legend()
plt.show()

plt.plot(temperatures,mui, label="mu_c(Ih)")
plt.plot(temperatures,muc0_g, label="mu_c0(sII)")
plt.plot(temperatures,mu_g, label="mu(sII)")
#plt.plot(temperatures,muw, label="mu_w")
plt.plot(temperatures,muw2, label="mu_w2")
plt.title("Chemical Potentials")
plt.xlabel("Temperature [K]")
plt.ylabel("F | U | H  [kJ/mol]")
plt.legend()
plt.show()

plt.plot(temperatures,hw, label="H(water)")
plt.plot(temperatures,hi, label="H(ice)")
plt.plot(temperatures,h, label="H(CH)")
plt.title("Enthalpy")
plt.xlabel("Temperature [K]")
plt.ylabel("F | U | H  [kJ/mol]")
plt.legend()
plt.show()

plt.plot(temperatures,hw-hi, label="\delta H(ice-water)")
plt.plot(temperatures,hw-h, label="\delta H(CH-water)")
plt.title("Latent heat")
plt.xlabel("Temperature [K]")
plt.ylabel("F | U | H  [kJ/mol]")
plt.legend()
plt.show()


