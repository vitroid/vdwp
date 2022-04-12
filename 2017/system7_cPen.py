#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#ふたたび通常のpythonに戻ってきた．
#ipython3 notebookは長いプログラムには向かない．

import histo2f
import histo
import derivative
import physconst as pc
import math
import errors
import numpy as np
import crystals
import sys
import interpolate as ip
import matplotlib.pyplot as plt
import os

def elimempty(values):
  newvalues = []
  for value in values:
    if value != "":
      newvalues.append(value)
  return newvalues


def usage():
  print("usage: {0} [-p][-T f,t,i][-L p,q,r][-G p,q,r][-X p,q,r][d]".format(sys.argv[0]))
  print("\t-p\t\tPressure [atm](1.0)")
  print("\t-T f,t,i\tTemperature range [From, To, Intv](263,303,3)")
  print("\t-L p,q,r\tLatent heat mode; p,q,r specifies molecular types in ID08.")
  print("\t-G p,q,r\tGraphics mode.")
  print("\t-X p,q,r\tTest mode.")
  print("\t-d\t\tDebug mode (False)")
  sys.exit(0)

    
    



#! Test case for methane hydrate
import molecule
def LoadMoleculeDict(filename):
    file = open(filename,encoding="utf-8")
    moldict = molecule.loadInfo( file )
    #Supplementary info
    for key,mol in moldict.items():
        if mol.moi[0] == 0:
            mol.symm = 1
            mol.dimen = 0
            mol.dof = 3
        elif mol.moi[2] == 0:
            mol.symm = 2      #not always correct
            mol.dimen = 1
            mol.dof = 5
        else:
            mol.symm = 2      #not always correct
            mol.dimen = 3
            mol.dof = 6
        mol.name = key
    #Aliases
    moldict["NEGTHF__"].name = "THF(invert)"
    moldict["CJTHF___"].name = "THF"
    moldict["CPENTANE"].name = "cPentane"
    moldict["CPENTAN+"].name = "cPentane (105%)"
    moldict["CPENTAN-"].name = "cPentane (95%)"
    moldict["CPEN+5L_"].name = "cPentane (105% L)"
    moldict["CPEN-5L_"].name = "cPentane (95% L)"
    moldict["CPENTA++"].name = "cPentane (110%)"
    moldict["CPENTA--"].name = "cPentane (90%)"
    moldict["LJME____"].name = "Methane"
    moldict["TYTHF___"].name = "THF(Yaga)"
    return moldict



moldict = LoadMoleculeDict("data/DEFR")

#
#変数の命名則
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

#User variables
pressure     = 101325.00 * 1# Pa
temperatures = np.array(np.arange(263,295,3))
mode = "Graphs"
mode = "None"
#guest = "TYTHF___"
#guest = "CPEN+5L_"
guest = "CPEN-5L_"
#guest = "CPENTANE"
#guest = "LJME____"
structure = "sII"
debug = False

if len(sys.argv) > 1:
    guest = sys.argv[1]


#分解状態での化学ポテンシャル(あるいは蒸気圧)
status = {"CJTHF___" : ("LIQUID-ANTOINE", 4.12118, 1202.942, -46.818),
          "NEGTHF__" : ("LIQUID-ANTOINE", 4.12118, 1202.942, -46.818),
          #"CPENTANE" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTANE" : ("LIQUID-LOGP"   ,   #進捗6より。
                        [ 252.1, 260., 270., 280., 290., 300. ],
                        [-2.49764352, -2.18736956, -1.72377413, -1.29552641, -0.80673168, -0.46915021]),
          "TIP4P   " : ("LIQUID-LOGP"   ,
                        [ 252.1, 260., 270., 280., 290., 300. ],
                        [-6.41141701, -5.95776395, -4.80618801, -4.00800805, -2.95803771, -2.85252098]),
          "CPENTAN+" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTAN-" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTA++" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTA--" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "LJME____" : ("IDEALGAS", pressure)}

def Raw2Table(lines,x,y):
    table = []
    for line in lines:
        cols = line.split()
        if len(cols) > max(x,y):
            table.append([float(cols[x]), float(cols[y])])
    return table


def RegisterStatus(label, FE, dens):
    ftemp,fes = zip(*FE)
    dtemp,dens = zip(*dens)
    if ftemp != dtemp:
        sys.exit(1)

    temp = np.array(dtemp)
    fes  = np.array(fes)
    dens = np.array(dens)
        
    #std pressure
    p0 = temp * 22.4 / 273.15
    #c-Pentane
    mol = moldict[label]

    pvap = p0*dens/mol.mass*np.exp(fes/(0.008314*temp))
    status[label] = ("LIQUID-LOGP", temp, np.log(pvap))
    

#/r2/matto/matto-gromacs/Pana/Neat/CPen+5E_/2015-01-04
Hneat = dict()
import os
dirs = os.listdir("MDdata/neat/1atm")
for g in dirs:
    if g[0] != ".":  #.DS_Store
        FE   = Raw2Table(open("MDdata/neat/1atm/"+g+"/FE.txt").readlines(), 0, 6)
        dens = Raw2Table(open("MDdata/neat/1atm/"+g+"/dens.txt").readlines(), 0, 1)
        RegisterStatus(g, FE, dens)
        #Enthalpy of the neat liquids
        H    = Raw2Table(open("MDdata/neat/1atm/"+g+"/H.txt").readlines(), 0, 1)
        temp,H = zip(*H)
        Hneat[g] = (temp,H)


#generate the list of guests, structures, and cages
cages = crystals.cage_components[structure][1]

import chempot

####### Clathrate-independent terms ####################################
#chemical potential of water ice
#Q構造のポテンシャルと振動の寄与から算出する。
import vdWP
errors.warn("Calculating chemical potential of water ice...")
mu_i = np.array([crystals.U_e["Ih"] + vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", T) for T in temperatures])
# h_i ###################### (50)
#h_i = mu_i - temperatures*derivative.derivatives(temperatures,mu_i)

#chemical potential of liquid water
#Antoineの式を使う。
import antoine

errors.warn("Calculating chemical potential of liquid water...")
#mu_w = np.array([ antoine.ChemicalPotentialOfWater(T) for T in temperatures ])
#運動エネルギーの項を含む場合。上とどうちがうのか。

water_mode = "ANTOINE" #or "GROMACS"
A = 5.40221
B = 1838.675
C = -31.737
tip4p = moldict["TIP4P   "]
TT    = status["TIP4P   "][1]
logP  = status["TIP4P   "][2]
if water_mode == "ANTOINE":
    #Antoine式のCPを使う
    mu_w = np.array([chempot.chempot(T,antoine.VaporPressure(T, A,B,C)) + chempot.IntegrationFix(T, tip4p.dimen)
             + chempot.StericFix(T, tip4p.mass, 2.0, tip4p.moi[0], tip4p.moi[1], tip4p.moi[2] )
             for T in temperatures])
elif water_mode == "GROMACS":
    #GROMACSによる純TIP4PのCPを使う
    mu_w = np.array([chempot.chempot(T,101326*math.exp(ip.value(T,TT,logP))) 
                         + chempot.IntegrationFix(T, tip4p.dimen)
                         + chempot.StericFix(T, tip4p.mass, 2.0, tip4p.moi[0], tip4p.moi[1], tip4p.moi[2] )
                         for T in temperatures])

    
if debug:
    plt.plot(temperatures,[math.exp(ip.value(T,TT,logP))*101326 for T in temperatures])
    plt.plot(temperatures,[antoine.VaporPressure(T,A,B,C) for T in temperatures])

#水のエンタルピー
#ここでは化学ポテンシャルの温度微分で計算している．
#GROMACSで直接計算することもできる．
#また，式52のように比熱を積分して求めることもできる．
# h_w ###################### (50)
h_w = mu_w - temperatures*derivative.derivatives(temperatures,mu_w)


#水と氷の273.15 KでのCPを求める．
mu_w_std = ip.value(273.15, temperatures, mu_w)
mu_i_std = ip.value(273.15, temperatures, mu_i)

adjust = "NONE"
#adjust = "DIFFERENTIAL"
#adjust = "PROPORTIONAL"
if adjust=="PROPORTIONAL":
    ratio = mu_w_std / mu_i_std
    #mu_i *= ratio
    mu_w /= ratio
elif adjust=="DIFFERENTIAL":
    diff = mu_w_std - mu_i_std
    #mu_i += diff
    mu_w -= diff
# h_i ###################### (50)
#ここで計算する．
h_i = mu_i - temperatures*derivative.derivatives(temperatures,mu_i)


#gas terms
mu_g = dict()

mol = moldict[guest]
errors.warn("Calculating chemical potential of guest "+mol.name+"...")
if guest not in status:
    errors.warn("Regard the guest in gas state at {0} Pa.".format(pressure))
    pressures = [pressure for T in temperatures]
else:
    if status[guest][0] == "IDEALGAS":
        errors.warn("Regard the guest in gas state at {0} Pa.".format(status[guest][1]))
        pressures = [status[guest][1] for T in temperatures]
    elif status[guest][0] == "LIQUID-ANTOINE":
        errors.warn("Regard the guest in ideal gas at its vapor pressure specified by Antoine's eq.")
        A = status[guest][1]
        B = status[guest][2]
        C = status[guest][3]
        pressures = [antoine.VaporPressure(T, A,B,C) for T in temperatures]
    elif status[guest][0] == "LIQUID-LOGP":
        errors.warn("Regard the guest in ideal gas at its vapor pressures given explicitly.")
        T    = status[guest][1]
        logP = status[guest][2]
        pressures = [101326*math.exp(ip.value(t, T, logP)) for t in temperatures]
# mu_g ######################(15)
mu_g = np.array([chempot.chempot(T,p) + chempot.IntegrationFix(T, mol.dimen) + chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2],debug=debug)  for T,p in zip(temperatures, pressures)])
# h_g ###################### (28)
h_g = dict()
mol = moldict[guest]
errors.warn("Calculating enthalpy of guest "+mol.name+"...")
h_g = mu_g - temperatures* np.array(derivative.derivatives(temperatures,mu_g))


####### Clathrate terms ################################################

#for hydrate structure types
errors.warn("Calculating chemical potential of empty clathrate "+structure+"...")
mu_e = np.array([ crystals.U_e[structure] + vdWP.FreeEnergyOfVibration(crystals.nma_file[structure], T) for T in temperatures ])
#水の自由エネルギーを操作する場合は，固体はいじらなくていい．
#if adjust == "PROPORTIONAL":
#    mu_e *= ratio
#if adjust == "DIFFERENTIAL":
#    mu_e += diff
    
#print mu_e
#空のケージの潜熱
# h_e ###################### (50)
errors.warn("Calculating enthalpy of empty clathrate "+structure+"...")
h_e = mu_e - temperatures* np.array(derivative.derivatives(temperatures,mu_e))


#for gases in the cage
f_c = dict()
h_c = dict()
for cage in cages:
    mol = moldict[guest]
    errors.warn("Calculating free energy of guest "+mol.name+" in a "+("%d" % cage)+" hedral cage...")

    histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
    histogram   = histo.loadAHisto( open(histofile) )

    # f_c ######################
    if histogram != None:
        f_c[cage] = np.array([histo2f.fvalue(histogram,T) + chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2]) for T in temperatures])
        # h_c ###################### (35)
        errors.warn("Calculating enthalpy of guest "+mol.name+" in a "+("%d" % cage)+" hedral cage...")
        h_c[cage] = f_c[cage] - temperatures* np.array(derivative.derivatives(temperatures,f_c[cage])) # +  (mol.dof/2.0)*pc.NkB*temperatures

#for gases in the lattice
mu_f = dict()
Deltamu = dict()
mol = moldict[guest]
errors.warn("Calculating chemical potential of guest "+mol.name+" in "+structure+"...")

Nw        = crystals.cage_components[structure][0]

# mu_f ######################
sum = np.zeros_like(temperatures)
beta = 1.0/(pc.NkB*temperatures)
for cage in f_c:
    coeffs = -pc.NkB*temperatures*cages[cage]/Nw
    body  = np.log(1+np.exp(beta*(mu_g-f_c[cage])))
    add = coeffs*body
    sum = sum + add
Deltamu = sum
mu_f = mu_e + sum
# h_f ###################### (49)
errors.warn("Calculating enthalpy of guest "+mol.name+" in "+structure+"...")
h_f = mu_f - temperatures* np.array(derivative.derivatives(temperatures,mu_f))


# If it is an ideal gas, calculation is much simpler in (31..33)
if debug:
  for guest in guests:
    #print "guest", guest
    #print "h_g",h_g[guest]
    dof = moldict[guest].dof
    x = (dof+2)/2.0
    #print "h_gi",temperatures * x *pc.NkB
    #print "ratio:",h_g[guest]/(temperatures * x *pc.NkB)
    #2015-2-4 confirmed for LJME____
    #2015-2-4 confirmed for LJPR2___
    #2015-2-4 confirmed for LJPRX___
    #No problem.

# Must be same as the averaged potential energy in the cage w_k in (36)
if debug:
  for guest,cage in guestandcage:
    histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
    histogram   = histo.loadAHisto( open(histofile) )
    if histogram != None:
      dof = moldict[guest].dof
      x = dof/2.0
      w_k = np.array([histo2f.energy(histogram,T) for T in temperatures]) + temperatures * x * pc.NkB
      #print "guest,cage", guest,cage
      #print "h_c",h_c[(guest,cage)]
      #print "w_k",w_k
      #print "ratio:",h_c[(guest,cage)]/w_k
      #print "diff:",h_c[(guest,cage)]- w_k
      #2015-2-5 confirmed for LJME____
      #2015-2-5 confirmed for LJPR2___
      #2015-2-5 confirmed for LJPRX___
      #No problem.





if debug:
    # h_2 ###################### (49)
    h_2 = np.zeros_like(temperatures)
    beta = 1.0/(pc.NkB*temperatures)
    Nw        = crystals.cage_components[structure][0]
    print("cages / Nw = {0s} / {1d}" % (cages,Nw))
    for cage in f_c:
        e = np.exp(beta*(mu_g - f_c[cage]))
        y = e/(1.0+e)
        h_2 = h_2 + y * cages[cage]/Nw* (h_g - h_c[cage])
    #Enthalpy of CH per water molecule
    h_f1 = h_e - h_2
    plt.plot(temperatures,h_f1)
    plt.plot(temperatures,h_f)
    #They are the same.

#Enthalpy of water and the free guest (51)
h_3 = np.zeros_like(temperatures)
for cage in f_c:
    e = np.exp(beta*(mu_g - f_c[cage]))
    y = e/(1.0+e)
    h_3 = h_3 + y * cages[cage]/Nw*h_g
# h_2 ###################### (51)
#田中先生のアドバイスによりh_3は加えない．
#h_d = h_w + h_3
h_d = h_w
#latent heat per guest molecule
#occupancy calculation
beta = 1.0/(pc.NkB*temperatures)
num_guest = np.zeros_like(temperatures)
for cage in f_c:
    e = np.exp(beta*(mu_g - f_c[cage]))
    y = e/(1.0+e)
    #errors.warn("occupancy in {0}:{1}".format(cage, y))
    num_guest = num_guest + y * cages[cage]/Nw
if num_guest[0] > 1e-4:
    if debug:
        print("Num of guest per a water molecule:", num_guest)
    mass_g = mol.mass
    mass_w = 18.0
    mass = mass_w + mol.mass*num_guest
    dh  = (h_d - h_f) / num_guest
    dhm = (h_d - h_f) / mass * 1000.0 #per kG
    dh_fwpg = dh      # delta H (filled - water) per guest
    dh_fwpm = dhm  # delta H (filled - water) per mass

#ここまでは解析解を主に利用した計算．










#以下ではMDの結果を積極的に利用する．
#まず，解析解とMDの結果を比較する．

#MDで求めた水のエンタルピー (GROMACS)
T_w = []
H_w = []
for line in """
252.1  -38.0179
260  -37.2569
270  -36.2945
280  -35.3782
290  -34.4608
300  -33.5675
""".split("\n"):
    cols = line.split()
    if 1 < len(cols):
        a,b = [float(x) for x in cols]
        T_w.append(a)
        H_w.append(b)
h_wMD = np.array([ip.value(T,T_w,H_w) for T in temperatures])
if mode == "Graphs":
    plt.figure(1)
    plt.plot(temperatures,h_w,     label="Water (Antoine)")
    plt.plot(temperatures,h_wMD,label="Gromacs TIP4P")
    plt.xlabel("Temperature [K]")
    plt.ylabel("H_w  [kJ/mol]")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))


#MDで求めたcPentaneのエンタルピー (GROMACS)
h_gMD = np.array([ip.value(T,*Hneat[guest]) for T in temperatures])
mol = moldict[guest]
if mode == "Graphs":
    plt.figure(2)
    plt.plot(temperatures,h_g,     label="{0} (from BAR)".format(mol.name))
    plt.plot(temperatures,h_gMD,label="Gromacs {0}".format(mol.name))
    plt.xlabel("Temperature [K]")
    plt.ylabel("H_g  [kJ/mol]")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))



cage=16
e = np.exp(beta*(mu_g - f_c[cage]))
y = e/(1.0+e)
num_guest = y * cages[cage]/Nw
#h_dMD = h_wMD + num_guest*h_gMD
#ここでも田中案によりゲスト項は加えない．
h_dMD = h_wMD
h_fMD = h_e - num_guest*(h_gMD - h_c[cage])

mol = moldict[guest]
mass_g = mol.mass
mass_w = 18.0
mass = mass_w + mass_g*num_guest
dh_fwpgMD  = (h_dMD - h_fMD) / num_guest
dh_fwpmMD  = (h_dMD - h_fMD) / mass * 1000.0 #per kG

plt.show()




#ここまでで一旦出力しようか．

#If you want to obtain the values at a specific temperature,
#you must set the temperature range at the beginning of this program in advance.
def TestMode(guest):
    mol = moldict[guest]
    T = temperatures[0]
    print("Test for",guest,"at",T,"K")
    beta = 1.0/(pc.NkB*T)
    print(mu_g[guest][0], "mu_g at ",pressure/101326,"atm")
    print(chempot.chempot(T,pressure), "chempot")
    print(chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2]),"StericFix")
    print("---")

    mol = moldict[guest]
    for cage in (12,14,16):
      histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
      histogram   = histo.loadAHisto( open(histofile) )
      if histogram:
        print(histo2f.fvalue(histogram,T), "raw f0_{0}".format(cage))
        print(chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2], debug=True), "StericFix for {0}hedral cage".format(cage))
        print("---")
    for cage in (14,12,16):
        if f_c.has_key((guest,cage)):
            print (f_c[(guest,cage)][0], "f_%d" % cage)
            print (mu_g[guest][0]-f_c[(guest,cage)][0], "mu_g - f_%d"%cage)
            print (math.exp(beta*(mu_g[guest][0]-f_c[(guest,cage)][0])), "exp(beta(mu_g - f_%d))"%cage)
    print (mu_e["sI"][0], "mu_c0 of empty clathrate sI")
    print (Deltamu[(guest,"sI")][0], "Delta mu_c sI")
    print (mu_f[(guest,"sI")][0], "mu_f sI")
    print (mu_e["sII"][0], "mu_c0 of empty clathrate sII")
    print (Deltamu[(guest,"sII")][0], "Delta mu_c sII")
    print (mu_f[(guest,"sII")][0], "mu_f sII")
    print()

def CageIntegrals(guest=""):
    if guest == "":
        print (("Cage Integrals "+"#"*72)[0:72])
        print ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format("guest",
                                                    "cage",
                                                    "f0_raw",
                                                    "ster",
                                                    "f",
                                                    "occup"))
    else:
        mol = moldict[guest]
        T = temperatures[0]
        for cage in (12,14,16):
            if cage in f_c:
                histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
                histogram   = histo.loadAHisto( open(histofile) )
                if histogram:
                    e = np.exp(beta*(mu_g-f_c[cage]))
                    y = e/(1.0+e)
                    print ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(guest, 
                                                                cage, 
                                                                histo2f.fvalue(histogram,T),
                                                                chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2]),
                                                                f_c[cage][0],
                                                                y[0]
                                                              ))


def Stabilities(guest=""):
    if guest == "":
        print(("Stabilities "+"#"*72)[0:72])
        print ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format("guest",
                                                         "struc",
                                                         "empty",
                                                         "guest",
                                                         "total",
                                                         "ice",
                                                         "rel_ice"))
    else:
        print ("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(guest,
                                                       structure,
                                                       mu_e[0],
                                                       Deltamu[0],
                                                       mu_f[0],
                                                       mu_i[0],
                                                       mu_f[0] - mu_i[0]))


def EmptyLattices(guest=""):
  print (("Empty Lattices "+"#"*72)[0:72])
  print ("{0}\t{1}\t{2}\t{3}".format("struc",
                                      "u",
                                      "f",
                                      "u+f"))
  T = temperatures[0]
  for struc in ("Ih","sI","sII"):
    if struc == "Ih":
      fvib = vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", T)      
    else:
      fvib = vdWP.FreeEnergyOfVibration(crystals.nma_file[struc], T)
    print ("{0}\t{1}\t{2}\t{3}".format(struc,
                                      crystals.U_e[struc],
                                      fvib,
                                      crystals.U_e[struc] + fvib))



def ChemicalPotentials(guest=""):
  if guest == "":
    print (("Chemical Potentials of the Guests "+"#"*72)[0:72])
    print ("{0}\t{1}".format("guest",
                            "mu_g"))
  else:
    print ("{0}\t{1}".format(guest,
                                                mu_g[0]))
  
    

def TestReport(guest):
    print (("Test Report "+"#"*72)[0:72])
    print ("Temperature:", temperatures[0])
    print ("Pressure:", pressure/101325)
    EmptyLattices()
    ChemicalPotentials()
    ChemicalPotentials(guest)
    CageIntegrals()
    CageIntegrals(guest)
    Stabilities()
    Stabilities(guest)
  

def GraphMode():
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['figure.figsize'] = (10.0, 8.0)
    font = {'family' : 'sans serif',
        'weight' : 'roman',
        'size'   : 14}
    matplotlib.rc('font', **font)
    plt.figure(1)
    plt.plot(temperatures,mu_i, label="mu_i ice Ih")
    plt.text(temperatures[0],mu_i[0], "ice Ih")
    plt.plot(temperatures,mu_e, label="empty " + structure)
    plt.text(temperatures[0],mu_e[0], "empty " + structure)
    #print(mu_e)
    #print(mu_f)
    mol = moldict[guest]
    plt.plot(temperatures,mu_f, label="mu_f "+structure+" "+mol.name)
    plt.text(temperatures[0],mu_f[0], structure+" "+mol.name)
    plt.plot(temperatures,mu_w, label="mu_w water")
    plt.text(temperatures[0],mu_w[0], "water")
    plt.title("Chemical Potential of Water in Various Phases")
    plt.xlabel("Temperature [K]")
    plt.ylabel("F  [kJ/mol]")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))

    plt.figure(2)
    plt.plot(temperatures,h_w, label="h_w water")
    plt.text(temperatures[0],h_w[0], "water")
    mol = moldict[guest]
    plt.plot(temperatures,h_d, label="h_d Sep. water+{0} ({1})".format(mol.name,structure))
    plt.text(temperatures[0],h_d[0], "Sep. water+{0} ({1})".format(mol.name,structure))
    plt.plot(temperatures,h_wMD, label="h_wMD water(MD)")
    plt.text(temperatures[0],h_wMD[0], "water(MD)")
    #plt.plot(temperatures,h_dMD, label="Sep. water+cPen(MD)")
    #plt.text(temperatures[0],h_dMD[0], "Sep. water+cPen(MD)")

    plt.plot(temperatures,h_fMD, label="h_fMD cPen CH(MD+vdWP)")
    plt.text(temperatures[0],h_fMD[0], "cPen CH(MD+vdWP)")



    plt.plot(temperatures,h_i, label="h_i ice Ih")
    plt.text(temperatures[0],h_i[0], "ice Ih")
    mol = moldict[guest]
    plt.plot(temperatures,h_f, label="h_f {0} {1} (vp+vdWP)".format(structure,mol.name))
    plt.text(temperatures[0],h_f[0], "{0} {1} (vp+vdWP)".format(structure,mol.name))
    plt.title("Enthalpy")
    plt.xlabel("Temperature [K]")
    plt.ylabel("H  [kJ/mol]")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))


    #per guest
    plt.figure(3)
    #for guest,structure in targets:
    mol = moldict[guest]
    plt.plot(temperatures,dh_fwpg, label=structure+" "+mol.name+"-water")
    plt.text(temperatures[0],dh_fwpg[0], structure+" "+mol.name+"-water")
    plt.plot(temperatures,dh_fwpgMD, label="sII cPen-water(MD)")
    plt.text(temperatures[0],dh_fwpgMD[0], "sII cPen-water(MD)")


    plt.title("Latent heat per guest molecule")
    plt.xlabel("Temperature [K]")
    plt.ylabel("\Delta H  [kJ/mol/guest mol]")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))

    #per kg
    plt.figure(4)
    plt.plot(temperatures,(h_w-h_i)/18.0*1000, label="ice-water")
    plt.text(temperatures[0],(h_w[0]-h_i[0])/18.0*1000, "ice-water")
    mol = moldict[guest]
    plt.plot(temperatures,dh_fwpm, label=structure+" "+mol.name+"-water")
    plt.text(temperatures[1],dh_fwpm[1], structure+" "+mol.name+"-water")
    plt.plot(temperatures,dh_fwpmMD, label="sII cPen-water(MD)")
    plt.text(temperatures[0],dh_fwpmMD[0], "sII cPen-water(MD)")


    plt.title("Latent heat per mass")
    plt.xlabel("Temperature [K]")
    plt.ylabel("\Delta H  [kJ/mol/kG]")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))

    plt.figure(5)
    e = np.exp(beta*(mu_g-f_c[16]))
    y = e/(1.0+e)
    plt.plot(temperatures, y, label="16-hedra")

    plt.title("Occupancy")
    plt.xlabel("Temperature [K]")
    plt.ylabel("Occupancy")
    plt.legend(loc='center left', bbox_to_anchor=(0.65,0.5))
    
    
    plt.show()


def LatentHeatMode():
    #only one guest
    for guest in guests:
        #determine melting point and crystal type
        mp = 0
        mp_struc = ""
        for structure in structures:
            mp_w = ip.zerocross(temperatures, mu_f[(guest, structure)] - mu_w3)
            if debug:
              print (mp_w,"MPW")
              print ("mu_f",mu_f[(guest,structure)])
              print ("mu_w",mu_w)
              print ("diff", mu_f[(guest,structure)] - mu_w)
            if temperatures[0] < mp_w < temperatures[-1]:
                if mp < mp_w:
                    mp_struc = structure
                    mp = mp_w
            # if temperatures[0] < mp_w < temperatures[-1]:
            #     if mp < mp_w and 273.15 < mp_w:
            #         mp_struc = structure
            #         mp = mp_w
            # mp_i = ip.zerocross(temperatures, mu_f[(guest, structure)] - mu_i)
            # if temperatures[0] < mp_i < temperatures[-1]:
            #     if mp < mp_i < 273.15:
            #         mp_struc = structure
            #         mp = mp_i
        if mp_struc != "":
          #latent = ip.zerocross( dh_fwpm[(guest, mp_struc)], temperatures - mp )
          latent = ip.value( mp, temperatures, dh_fwpm )
          print("{0}\t{1}\t{2}\t{3}\t# kJ/kG".format(guest,mp_struc, mp,latent))

if mode == "Test":
  TestReport(guest)

if mode == "Graphs":
  GraphMode()

if mode == "Latent":
  LatentHeatMode()
print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format("#temp", "h_gMD", "mu_g", "f_16", "mu-f", "mu_f", "name"))
print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(270, ip.value( 270, temperatures, h_gMD ), ip.value( 270, temperatures, mu_g ), ip.value( 270, temperatures, f_c[16] ), ip.value( 270, temperatures, mu_g - f_c[16] ), ip.value( 270, temperatures, mu_f ), mol.name ))

