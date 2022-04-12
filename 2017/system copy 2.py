#!/usr/bin/env python
#coding: utf-8

import histo2f
import histo
import derivative
import physconst as pc
import math
import errors

#とりあえず、pretty printingはあとまわしにして、すべての計算をpython上で行えるようにする。

#! Test case for methane hydrate
import molecule
file = open("data/DEFR")
moldict = molecule.loadInfo( file )

#
#f, mu: free energy/chemical potential
#h    : enthalpy
#
#_w   : of water
#_i   : of ice
#_e   : of the empty lattice 
#_f   : of the filled lattice 
#_g   : of the guest
#_c   : in the cage
#


pressure     = 101325.00 # Pa
structures = {"sII" : (136,{16:8.,12:16.}), #cage 0 is Nw
              "sI"  : (46, {12:2.,14:6.})}
nma_file = {"sII" : "data/C15retry/cat.0.renma.dist",
            "sI"  : "data/A15retry/cat.0.renma.dist" }
U_e      = {"sII" :  -54.7751,
            "sI"  :  -54.3857}
status = {"CJTHF___" : ("LIQUID-ANTOINE",4.12118, 1202.942,-46.818),
          "CPENTANE" : ("LIQUID-ANTOINE",4.24714, 1235.305, -30.666),
          "CPENTAN+" : ("LIQUID-ANTOINE",4.24714, 1235.305, -30.666),
          "CPENTAN-" : ("LIQUID-ANTOINE",4.24714, 1235.305, -30.666),
          "CPENTA++" : ("LIQUID-ANTOINE",4.24714, 1235.305, -30.666),
          "CPENTA--" : ("LIQUID-ANTOINE",4.24714, 1235.305, -30.666),
          "LJME____" : ("IDEALGAS", pressure)}
#Supplementary info
moldict["CJTHF___"].symm = 2
moldict["CPENTANE"].symm = 2
moldict["CPENTA++"].symm = 2
moldict["CPENTA--"].symm = 2
moldict["CPENTAN+"].symm = 2
moldict["CPENTAN-"].symm = 2
moldict["LJME____"].symm = 1
moldict["CJTHF___"].dimen = 3
moldict["CPENTANE"].dimen = 3
moldict["CPENTAN+"].dimen = 3
moldict["CPENTAN-"].dimen = 3
moldict["CPENTA++"].dimen = 3
moldict["CPENTA--"].dimen = 3
moldict["LJME____"].dimen = 0
moldict["CJTHF___"].dof = 6
moldict["CPENTANE"].dof = 6
moldict["CPENTAN+"].dof = 6
moldict["CPENTAN-"].dof = 6
moldict["CPENTA++"].dof = 6
moldict["CPENTA--"].dof = 6
moldict["LJME____"].dof = 3
moldict["CJTHF___"].name = "THF"
moldict["CPENTANE"].name = "cPentane"
moldict["CPENTAN+"].name = "cPentane (105%)"
moldict["CPENTAN-"].name = "cPentane (95%)"
moldict["CPENTA++"].name = "cPentane (110%)"
moldict["CPENTA--"].name = "cPentane (90%)"
moldict["LJME____"].name = "Methane"


import numpy

temperatures = numpy.array(range(273-30,283,3))
#guests = ("CJTHF___", "CPENTANE", "CPENTAN+", "CPENTAN-")
#guests = ("CPENTA++", "CPENTAN+", "CPENTANE", "CPENTAN-","CPENTA--")
targets = (("CPENTA++","sII"),)



import chempot

#chemical potential of water ice
#Q構造のポテンシャルと振動の寄与から算出する。
import vdWP
errors.warn("Calculating chemical potentials of water ice...")
mu_i = [ -55.2238 + vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", T) for T in temperatures ]
h_i = mu_i - temperatures*derivative.derivatives(temperatures,mu_i)

#chemical potential of liquid water
#Antoineの式を使う。
import antoine

errors.warn("Calculating chemical potentials of liquid water...")
mu_w = [ antoine.ChemicalPotentialOfWater(T) for T in temperatures ]
#運動エネルギーの項を含む場合。上とどうちがうのか。
A = 5.40221
B = 1838.675
C = -31.737
tip4p = moldict["TIP4P   "]
mu_w2 = [chempot.RigidRotor(T,antoine.VaporPressure(T, A,B,C), tip4p.moi[0], tip4p.moi[1], tip4p.moi[2], tip4p.mass, 2.0) for T in temperatures]
h_w = mu_w2 - temperatures*derivative.derivatives(temperatures,mu_w2)


#for hydrate structures
mu_e = dict()
for structure in U_e:
    errors.warn("Calculating chemical potentials of empty clathrate "+structure+"...")
    mu_e[structure] = [ U_e[structure] + vdWP.FreeEnergyOfVibration(nma_file[structure], T) for T in temperatures ]
#print mu_e

#for Guests
mu_g = dict()
mu_f = dict()
f_c = dict()
Deltamu = dict()
for guest in guests:
    mol = moldict[guest]
    errors.warn("Calculating chemical potentials of guest "+mol.name+"...")
    if status[guest][0] == "IDEALGAS":
        pressures = [status[guest][1] for T in temperatures]
    if status[guest][0] == "LIQUID-ANTOINE":
        A = status[guest][1]
        B = status[guest][2]
        C = status[guest][3]
        pressures = [antoine.VaporPressure(T, A,B,C) for T in temperatures]

    if mol.dimen == 3:
        mu_g[guest] = numpy.array([chempot.RigidRotor(temperatures[i], pressures[i], mol.moi[0], mol.moi[1], mol.moi[2], mol.mass, mol.symm) for i in range(len(temperatures))])
    elif mol.dimen == 1:
        errors.warn("Rodlike molecule is not implemented for now.")
    elif mol.dimen == 0:
        mu_g[guest] = numpy.array([chempot.Monatomic(temperatures[i], pressures[i], mol.mass) for i in range(len(temperatures))])


    structure = structures[guest][0]
    Nw        = structures[guest][1]
    cages     = structures[guest][2]

    histogram = dict()
    for cage in cages:
        histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
        histogram[cage] = histo.loadAHisto( open(histofile) )

    for cage in cages:
        if mol.dimen == 3:
            f_c[(guest,cage)] = numpy.array([histo2f.fvalue(histogram[cage],T)+chempot.Avalue(T,mol.mass)+chempot.Cvalue(T,mol.moi[0], mol.moi[1], mol.moi[2])+chempot.SymmetryFix(T,mol.symm)+chempot.IntegrationFix(T,mol.dimen) for T in temperatures])
        elif mol.dimen == 1: #??
            errors.warn("Rodlike molecule is not implemented for now.")
        elif mol.dimen == 0:
            f_c[(guest,cage)] = numpy.array([histo2f.fvalue(histogram[cage],T)+chempot.Avalue(T,mol.mass)+chempot.SymmetryFix(T,mol.symm)+chempot.IntegrationFix(T,mol.dimen) for T in temperatures])


    Deltamu[guest] = []
    for i in range(len(temperatures)):
        T = temperatures[i]
        beta = 1.0/(pc.NkB*T)
        sum = 0.0
        for cage in cages:
            sum += (-pc.NkB*T*cages[cage]/Nw*math.log(1+math.exp(beta*(mu_g[guest][i]-f_c[(guest,cage)][i]))))
        Deltamu[guest].append(sum)
    Deltamu[guest] = numpy.array(Deltamu[guest])
    mu_f[guest] = mu_e[structure] + Deltamu[guest]

#結局、実験の蒸気圧をつかった数字がのきなみ信用できない。
#それはさておき、潜熱の計算にいよいよとりかかる。

#空のケージの潜熱
h_e = dict()
for structure in U_e:
    errors.warn("Calculating enthalpies of empty clathrate "+structure+"...")
    h_e[structure] = mu_e[structure] - temperatures* numpy.array(derivative.derivatives(temperatures,mu_e[structure]))


h_g = dict()
h_c = dict()
h_f = dict()

dh_fwpg = dict() # delta H (filled - water) per guest
dh_fwpm = dict() # delta H (filled - water) per mass
for guest in guests:
    mol = moldict[guest]
    errors.warn("Calculating enthalpies of guest "+mol.name+"...")
    h_g[guest] = mu_g[guest] - temperatures* numpy.array(derivative.derivatives(temperatures,mu_g[guest]))
    structure = structures[guest][0]
    Nw        = structures[guest][1]
    cages     = structures[guest][2]

    for cage in cages:
        h_c[(guest,cage)] = f_c[(guest,cage)] - temperatures* numpy.array(derivative.derivatives(temperatures,f_c[(guest,cage)])) +  (mol.dof/2.0)*pc.NkB*temperatures


    DeltaH = []
    for i in range(len(temperatures)):
        T = temperatures[i]
        beta = 1.0/(pc.NkB*T)
        sum = 0.0
        for cage in cages:
            e = math.exp(beta*(mu_g[guest][i]-f_c[(guest,cage)][i]))
            y = e/(1.0+e)
            sum += y * cages[cage]/Nw* (h_g[guest][i] - h_c[(guest,cage)][i])
        DeltaH.append(sum)
    DeltaH = numpy.array(DeltaH)
    #latent heat per water molecule
    h_f[guest] = h_e[structure] - DeltaH

    #latent heat per guest molecule
    dh = []
    dhm = []
    for i in range(len(temperatures)):
        T = temperatures[i]
        beta = 1.0/(pc.NkB*T)
        sum = 0.0
        for cage in cages:
            e = math.exp(beta*(mu_g[guest][i]-f_c[(guest,cage)][i]))
            y = e/(1.0+e)
            sum += y * cages[cage]/Nw
        mass_g = mol.mass
        mass_w = 18.0
        mass = mass_w + mol.mass*sum
        dh.append( (h_w[i] - h_f[guest][i]) / sum )
        dhm.append( (h_w[i] - h_f[guest][i]) / mass * 1000.0 ) #per kG
    dh_fwpg[guest] = numpy.array(dh)
    dh_fwpm[guest] = numpy.array(dhm)


def test(guest):
    T = temperatures[0]
    beta = 1.0/(pc.NkB*T)
    print f_c[(guest,14)][0], "Delta f_14"
    print math.exp(beta*(mu_g[guest][0]-f_c[(guest,14)][0])), "exp(beta(mu-f14))"
    print f_c[(guest,12)][0], "Delta f_12"
    print math.exp(beta*(mu_g[guest][0]-f_c[(guest,12)][0])), "exp(beta(mu-f12))"
    print Deltamu[guest][0], "Delta mu_c"
    print mu_f[guest][0], "mu_f"

#test("LJME____")



import matplotlib.pyplot as plt

plt.figure(1)
plt.plot(temperatures,mu_i, label="ice Ih")
for structure in U_e:
    plt.plot(temperatures,mu_e[structure], label="empty " + structure)
for guest in guests:
    structure = structures[guest][0]
    mol = moldict[guest]
    plt.plot(temperatures,mu_f[guest], label=structure+" "+mol.name)
plt.plot(temperatures,mu_w2, label="water")
plt.title("Chemical Potentials")
plt.xlabel("Temperature [K]")
plt.ylabel("F  [kJ/mol]")
plt.legend()

plt.figure(2)
plt.plot(temperatures,h_w, label="water")
plt.plot(temperatures,h_i, label="ice Ih")
for guest in guests:
    structure = structures[guest][0]
    mol = moldict[guest]
    plt.plot(temperatures,h_f[guest], label=structure+" "+mol.name)
plt.title("Enthalpy")
plt.xlabel("Temperature [K]")
plt.ylabel("H  [kJ/mol]")
plt.legend()

plt.figure(3)
plt.plot(temperatures,h_w-h_i, label="ice-water")
for guest in guests:
    structure = structures[guest][0]
    mol = moldict[guest]
    plt.plot(temperatures,h_w-h_f[guest], label=structure+" "+mol.name+"-water")
#    plt.plot(temperatures,h_i-h_f[guest], label=structure+" "+mol.name+"-ice")
plt.title("Latent heat per water molecule")
plt.xlabel("Temperature [K]")
plt.ylabel("\Delta H  [kJ/mol/water mol]")
plt.legend()

#per guest
plt.figure(4)
for guest in guests:
    structure = structures[guest][0]
    mol = moldict[guest]
    plt.plot(temperatures,dh_fwpg[guest], label=structure+" "+mol.name+"-water")
plt.title("Latent heat per guest molecule")
plt.xlabel("Temperature [K]")
plt.ylabel("\Delta H  [kJ/mol/guest mol]")
plt.legend()

#per kg
plt.figure(5)
plt.plot(temperatures,(h_w-h_i)/18.0*1000, label="ice-water")
for guest in guests:
    structure = structures[guest][0]
    mol = moldict[guest]
    plt.plot(temperatures,dh_fwpm[guest], label=structure+" "+mol.name+"-water")
plt.title("Latent heat per mass")
plt.xlabel("Temperature [K]")
plt.ylabel("\Delta H  [kJ/mol/kG]")
plt.legend()

plt.show()

