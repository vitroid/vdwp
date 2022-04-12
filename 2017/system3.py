#!/usr/bin/env python
#coding: utf-8

import histo2f
import histo
import derivative
import physconst as pc
import math
import errors
import numpy
import crystals
import sys

#とりあえず、pretty printingはあとまわしにして、すべての計算をpython上で行えるようにする。
#時間のかかる計算をnumpyの関数に置きかえ、時間を測定する。(ヒストグラムの積分?)

#2014-12 田中さんとの議論では、剛体分子だけ結果がおかしいなら座標変換の問題ではないか、ということ。
#ただ、回転の順番を間違えたとしても、全方向をサンプルしていれば問題ないような気もする。
#とりあえず、yaplotで分子配向を表示してみる。

#floating point range iterator
def frange(x, y, jump):
  while x < y:
    yield x
    x += jump

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

#User variables
pressure     = 101325.00 * 1# Pa
temperatures = numpy.array(range(273-10,273+30,3))
mode = "Graphs"
debug = False

def usage():
  print "usage: {} [-p][-T f,t,i][-L p,q,r][-G p,q,r][-X p,q,r][d]".format(sys.argv[0])
  print "\t-p\t\tPressure [atm](1.0)"
  print "\t-T f,t,i\tTemperature range [From, To, Intv](263,303,3)"
  print "\t-L p,q,r\tLatent heat mode; p,q,r specifies molecular types in ID08."
  print "\t-G p,q,r\tGraphics mode."
  print "\t-X p,q,r\tTest mode."
  print "\t-d\t\tDebug mode (False)"
  sys.exit(0)


def elimempty(values):
  newvalues = []
  for value in values:
    if value != "":
      newvalues.append(value)
  return newvalues



while len(sys.argv) > 1:
    arg = sys.argv.pop(1)
    if arg == "-h":
      usage()
    elif arg == "-p":
        pressure = float(sys.argv.pop(1)) * 101326.0
    elif arg == "-d":
        debug = True
    elif arg == "-T":
        values = sys.argv.pop(1).split(",")
        values = [float(x) for x in values]
        temperatures = numpy.array(list(frange(*values)))
    elif arg == "-L":  #for those interested only on latent heats
        mode = "Latent"
        guests = elimempty(sys.argv.pop(1).split(","))
        targets = [(guest,structure) for structure in crystals.cage_components.keys() for guest in guests]
    elif arg == "-G":  #graph
        mode = "Graphs"
        guests = elimempty(sys.argv.pop(1).split(","))
        targets = [(guest,structure) for structure in crystals.cage_components.keys() for guest in guests]
    elif arg == "-X":  #test
        mode = "Test"
        guests = elimempty(sys.argv.pop(1).split(","))
        targets = [(guest,structure) for structure in crystals.cage_components.keys() for guest in guests]
    else:
        print("Ignored argument: %s".format(arg))


status = {"CJTHF___" : ("LIQUID-ANTOINE", 4.12118, 1202.942, -46.818),
          "NEGTHF__" : ("LIQUID-ANTOINE", 4.12118, 1202.942, -46.818),
          "CPENTANE" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTAN+" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTAN-" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTA++" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "CPENTA--" : ("LIQUID-ANTOINE", 4.24714, 1235.305, -30.666),
          "LJME____" : ("IDEALGAS", pressure)}
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
moldict["CPENTA++"].name = "cPentane (110%)"
moldict["CPENTA--"].name = "cPentane (90%)"
moldict["LJME____"].name = "Methane"

#generate the list of guests, structures, and cages
guestandcage = set()
structures = set()
for guest,structure in targets:
    structures.add(structure)
    for cage in crystals.cage_components[structure][1]:
        guestandcage.add((guest,cage))

import chempot

####### Clathrate-independent terms ####################################
#chemical potential of water ice
#Q構造のポテンシャルと振動の寄与から算出する。
import vdWP
errors.warn("Calculating chemical potential of water ice...")
mu_i = numpy.array([ crystals.U_e["Ih"] + vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", T) for T in temperatures ])
h_i = mu_i - temperatures*derivative.derivatives(temperatures,mu_i)

#chemical potential of liquid water
#Antoineの式を使う。
import antoine

errors.warn("Calculating chemical potential of liquid water...")
mu_w = numpy.array([ antoine.ChemicalPotentialOfWater(T) for T in temperatures ])
#運動エネルギーの項を含む場合。上とどうちがうのか。
A = 5.40221
B = 1838.675
C = -31.737
tip4p = moldict["TIP4P   "]
mu_w2 = numpy.array([chempot.chempot(T,antoine.VaporPressure(T, A,B,C)) + chempot.IntegrationFix(T, tip4p.dimen)
         + chempot.StericFix(T, tip4p.mass, 2.0, tip4p.moi[0], tip4p.moi[1], tip4p.moi[2] )
         for T in temperatures])
h_w = mu_w2 - temperatures*derivative.derivatives(temperatures,mu_w2)

offset = ( crystals.U_e["Ih"] + vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", 273)
           - (chempot.chempot(273,antoine.VaporPressure(273, A,B,C)) + chempot.IntegrationFix(273, tip4p.dimen)
         + chempot.StericFix(273, tip4p.mass, 2.0, tip4p.moi[0], tip4p.moi[1], tip4p.moi[2] )))
mu_w3 = mu_w2 + offset
#mu_w3 = mu_w2
h_w = mu_w3 - temperatures*derivative.derivatives(temperatures,mu_w3)

#THIS IS TEST.
#mu_i -= offset
#h_i = mu_i - temperatures*derivative.derivatives(temperatures,mu_i)



#gas terms
mu_g = dict()

for guest in guests:
    mol = moldict[guest]
    errors.warn("Calculating chemical potential of guest "+mol.name+"...")
    if not status.has_key(guest):
        pressures = [pressure for T in temperatures]
    else:
        if status[guest][0] == "IDEALGAS":
            pressures = [status[guest][1] for T in temperatures]
        if status[guest][0] == "LIQUID-ANTOINE":
            A = status[guest][1]
            B = status[guest][2]
            C = status[guest][3]
            pressures = [antoine.VaporPressure(T, A,B,C) for T in temperatures]

    # mu_g ######################
    mu_g[guest] = numpy.array([chempot.chempot(T,p) + chempot.IntegrationFix(T, mol.dimen) + chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2],debug=debug)  for T,p in zip(temperatures, pressures)])


####### Clathrate terms ################################################

#for hydrate structure types
mu_e = dict()
for structure in structures:
    errors.warn("Calculating chemical potential of empty clathrate "+structure+"...")
    mu_e[structure] = numpy.array([ crystals.U_e[structure] + vdWP.FreeEnergyOfVibration(crystals.nma_file[structure], T) for T in temperatures ])
#print mu_e


#for gases in the cage
f_c = dict()
for guest,cage in guestandcage:
    mol = moldict[guest]
    errors.warn("Calculating free energy of guest "+mol.name+" in a "+("%d" % cage)+" hedral cage...")

    histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
    histogram   = histo.loadAHisto( open(histofile) )

    # f_c ######################
    if histogram != None:
        f_c[(guest,cage)] = numpy.array([histo2f.fvalue(histogram,T) + chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2]) for T in temperatures])

#for gases in the lattice
mu_f = dict()
Deltamu = dict()
for guest,structure in targets:
    mol = moldict[guest]
    errors.warn("Calculating chemical potential of guest "+mol.name+" in "+structure+"...")

    Nw        = crystals.cage_components[structure][0]
    cages     = crystals.cage_components[structure][1]

    # mu_f ######################
    Deltamu[(guest,structure)] = []
    for i in range(len(temperatures)):
        T = temperatures[i]
        beta = 1.0/(pc.NkB*T)
        sum = 0.0
        for cage in cages:
            if f_c.has_key((guest,cage)):
                sum += (-pc.NkB*T*cages[cage]/Nw*math.log(1+math.exp(beta*(mu_g[guest][i]-f_c[(guest,cage)][i]))))
        Deltamu[(guest,structure)].append(sum)
    Deltamu[(guest,structure)] = numpy.array(Deltamu[(guest,structure)])
    mu_f[(guest,structure)] = mu_e[structure] + Deltamu[(guest,structure)]

#結局、実験の蒸気圧をつかった数字がのきなみ信用できない。
#それはさておき、潜熱の計算にいよいよとりかかる。

#空のケージの潜熱
# h_e ###################### (50)
h_e = dict()
for structure in structures:
    errors.warn("Calculating enthalpy of empty clathrate "+structure+"...")
    h_e[structure] = mu_e[structure] - temperatures* numpy.array(derivative.derivatives(temperatures,mu_e[structure]))

# h_g ###################### (28)
h_g = dict()
for guest in guests:
    mol = moldict[guest]
    errors.warn("Calculating enthalpy of guest "+mol.name+"...")
    h_g[guest] = mu_g[guest] - temperatures* numpy.array(derivative.derivatives(temperatures,mu_g[guest]))
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

# h_c ###################### (35)
h_c = dict()
#for gases in the cage
for guest,cage in guestandcage:
    mol = moldict[guest]
    errors.warn("Calculating enthalpy of guest "+mol.name+" in a "+("%d" % cage)+" hedral cage...")
    if f_c.has_key((guest,cage)):
        h_c[(guest,cage)] = f_c[(guest,cage)] - temperatures* numpy.array(derivative.derivatives(temperatures,f_c[(guest,cage)])) # +  (mol.dof/2.0)*pc.NkB*temperatures
# Must be same as the averaged potential energy in the cage w_k in (36)
if debug:
  for guest,cage in guestandcage:
    histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
    histogram   = histo.loadAHisto( open(histofile) )
    if histogram != None:
      dof = moldict[guest].dof
      x = dof/2.0
      w_k = numpy.array([histo2f.energy(histogram,T) for T in temperatures]) + temperatures * x * pc.NkB
      #print "guest,cage", guest,cage
      #print "h_c",h_c[(guest,cage)]
      #print "w_k",w_k
      #print "ratio:",h_c[(guest,cage)]/w_k
      #print "diff:",h_c[(guest,cage)]- w_k
      #2015-2-5 confirmed for LJME____
      #2015-2-5 confirmed for LJPR2___
      #2015-2-5 confirmed for LJPRX___
      #No problem.


# h_f in (49).
h_f = dict()

dh_fwpg = dict() # delta H (filled - water) per guest
dh_fwpm = dict() # delta H (filled - water) per mass
for guest,structure in targets:
    mol = moldict[guest]
    errors.warn("Calculating enthalpy of guest "+mol.name+" in "+structure+"...")

    Nw        = crystals.cage_components[structure][0]
    cages     = crystals.cage_components[structure][1]
    if debug:
      print "cages / Nw = %s / %d" % (cages,Nw)
    # h_2 ###################### (49)
    h_2 = numpy.zeros(len(temperatures))
    beta = 1.0/(pc.NkB*temperatures)
    for cage in cages:
      if f_c.has_key((guest,cage)):
        e = numpy.exp(beta*(mu_g[guest]-f_c[(guest,cage)]))
        y = e/(1.0+e)
        h_2 += y * cages[cage]/Nw* (h_g[guest] - h_c[(guest,cage)])
    #Enthalpy of CH per water molecule
    h_f[(guest,structure)] = h_e[structure] - h_2
    #Enthalpy of water and the free guest (51)
    h_3 = numpy.zeros(len(temperatures))
    for cage in cages:
      if f_c.has_key((guest,cage)):
        e = numpy.exp(beta*(mu_g[guest]-f_c[(guest,cage)]))
        y = e/(1.0+e)
        h_3 += y * cages[cage]/Nw*h_g[guest]
    h_d = h_w + h_3
    #latent heat per guest molecule
    #occupancy calculation
    beta = 1.0/(pc.NkB*temperatures)
    num_guest = numpy.zeros(len(temperatures))
    for cage in cages:
      if f_c.has_key((guest,cage)):
        e = numpy.exp(beta*(mu_g[guest]-f_c[(guest,cage)]))
        y = e/(1.0+e)
        if debug:
          print "occupancy in %d:%s" % (cage, y)
        num_guest += y * cages[cage]/Nw
    if num_guest[0] != 0.0:
      if debug:
        print "Num of guest per a water molecule:", num_guest
      mass_g = mol.mass
      mass_w = 18.0
      mass = mass_w + mol.mass*num_guest
      dh  = (h_d - h_f[(guest,structure)]) / num_guest
      dhm = (h_d - h_f[(guest,structure)]) / mass * 1000.0 #per kG
      dh_fwpg[(guest,structure)] = dh
      dh_fwpm[(guest,structure)] = dhm
      if debug:
        print "dhm:", dhm
        print "LH of ice:",  (h_w - h_i) / mass_w * 1000


#If you want to obtain the values at a specific temperature,
#you must set the temperature range at the beginning of this program in advance.
def TestMode(guest):
    mol = moldict[guest]
    T = temperatures[0]
    print "Test for",guest,"at",T,"K"
    beta = 1.0/(pc.NkB*T)
    print mu_g[guest][0], "mu_g at ",pressure/101326,"atm"
    print chempot.chempot(T,pressure), "chempot"
    print chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2]),"StericFix"
    print "---"

    mol = moldict[guest]
    for cage in (12,14,16):
      histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
      histogram   = histo.loadAHisto( open(histofile) )
      if histogram:
        print histo2f.fvalue(histogram,T), "raw f0_{0}".format(cage)
        print chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2], debug=True), "StericFix for {0}hedral cage".format(cage)
        print "---"
    for cage in (14,12,16):
        if f_c.has_key((guest,cage)):
            print f_c[(guest,cage)][0], "f_%d" % cage
            print mu_g[guest][0]-f_c[(guest,cage)][0], "mu_g - f_%d"%cage
            print math.exp(beta*(mu_g[guest][0]-f_c[(guest,cage)][0])), "exp(beta(mu_g - f_%d))"%cage
    print mu_e["sI"][0], "mu_c0 of empty clathrate sI"
    print Deltamu[(guest,"sI")][0], "Delta mu_c sI"
    print mu_f[(guest,"sI")][0], "mu_f sI"
    print mu_e["sII"][0], "mu_c0 of empty clathrate sII"
    print Deltamu[(guest,"sII")][0], "Delta mu_c sII"
    print mu_f[(guest,"sII")][0], "mu_f sII"
    print

def CageIntegrals(guest=""):
  if guest == "":
    print ("Cage Integrals "+"#"*72)[0:72]
    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format("guest",
                                                "cage",
                                                "f0_raw",
                                                "ster",
                                                "f",
                                                "occup")
  else:
    mol = moldict[guest]
    T = temperatures[0]
    for cage in (12,14,16):
      histofile   = "data/" + guest + "." + ("%d" % cage) + "hedra.histo"
      histogram   = histo.loadAHisto( open(histofile) )
      if histogram:
        e = numpy.exp(beta*(mu_g[guest]-f_c[(guest,cage)]))
        y = e/(1.0+e)
        print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(guest, 
                                                    cage, 
                                                    histo2f.fvalue(histogram,T),
                                                    chempot.StericFix(T, mol.mass, mol.symm, mol.moi[0], mol.moi[1], mol.moi[2]),
                                                    f_c[(guest,cage)][0],
                                                    y[0]
                                                  )


def Stabilities(guest=""):
  if guest == "":
    print ("Stabilities "+"#"*72)[0:72]
    print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format("guest",
                                                     "struc",
                                                     "empty",
                                                     "guest",
                                                     "total",
                                                     "ice",
                                                     "rel_ice")
  else:
    for struc in sorted(("sI","sII"), cmp=lambda x,y: cmp(mu_f[(guest,x)][0],mu_f[(guest,y)][0])):
      print "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}".format(guest,
                                                       struc,
                                                       mu_e[struc][0],
                                                       Deltamu[(guest,struc)][0],
                                                       mu_f[(guest,struc)][0],
                                                       mu_i[0],
                                                       mu_f[(guest,struc)][0] - mu_i[0])


def EmptyLattices(guest=""):
  print ("Empty Lattices "+"#"*72)[0:72]
  print "{0}\t{1}\t{2}\t{3}".format("struc",
                                      "u",
                                      "f",
                                      "u+f")
  T = temperatures[0]
  for struc in ("Ih","sI","sII"):
    if struc == "Ih":
      fvib = vdWP.FreeEnergyOfVibration("data/Ih-Tanaka1998retry/cat.0.renma.dist", T)      
    else:
      fvib = vdWP.FreeEnergyOfVibration(crystals.nma_file[struc], T)
    print "{0}\t{1}\t{2}\t{3}".format(struc,
                                      crystals.U_e[struc],
                                      fvib,
                                      crystals.U_e[struc] + fvib)



def ChemicalPotentials(guest=""):
  if guest == "":
    print ("Chemical Potentials of the Guests "+"#"*72)[0:72]
    print "{0}\t{1}".format("guest",
                            "mu_g")
  else:
    print "{0}\t{1}".format(guest,
                                                mu_g[guest][0])
  
    

def TestReport(guests):
  print ("Test Report "+"#"*72)[0:72]
  print "Temperature:", temperatures[0]
  print "Pressure:", pressure/101325
  EmptyLattices()
  ChemicalPotentials()
  for guest in guests:
    ChemicalPotentials(guest)
  CageIntegrals()
  for guest in guests:
    CageIntegrals(guest)
  Stabilities()
  for guest in guests:
    Stabilities(guest)
  

def GraphMode():
    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.plot(temperatures,mu_i, label="ice Ih")
    plt.text(temperatures[0],mu_i[0], "ice Ih")
    for structure in structures:
        plt.plot(temperatures,mu_e[structure], label="empty " + structure)
        plt.text(temperatures[0],mu_e[structure][0], "empty " + structure)
    for guest,structure in targets:
        mol = moldict[guest]
        plt.plot(temperatures,mu_f[(guest,structure)], label=structure+" "+mol.name)
        plt.text(temperatures[0],mu_f[(guest,structure)][0], structure+" "+mol.name)
    plt.plot(temperatures,mu_w3, label="water")
    plt.text(temperatures[0],mu_w3[0], "water")
    plt.title("Chemical Potential")
    plt.xlabel("Temperature [K]")
    plt.ylabel("F  [kJ/mol]")
    plt.legend()

    plt.figure(2)
    plt.plot(temperatures,h_w, label="water")
    plt.text(temperatures[0],h_w[0], "water")
    plt.plot(temperatures,h_i, label="ice Ih")
    plt.text(temperatures[0],h_i[0], "ice Ih")
    for guest,structure in targets:
        mol = moldict[guest]
        plt.plot(temperatures,h_f[(guest,structure)], label=structure+" "+mol.name)
        plt.text(temperatures[0],h_f[(guest,structure)][0], structure+" "+mol.name)
    plt.title("Enthalpy")
    plt.xlabel("Temperature [K]")
    plt.ylabel("H  [kJ/mol]")
    plt.legend()

    plt.figure(3)
    plt.plot(temperatures,h_w-h_i, label="ice-water")
    plt.text(temperatures[0],h_w[0]-h_i[0], "ice-water")
    for guest,structure in targets:
        mol = moldict[guest]
        plt.plot(temperatures,h_w-h_f[(guest,structure)], label=structure+" "+mol.name+"-water")
        plt.text(temperatures[0],h_w[0]-h_f[(guest,structure)][0], structure+" "+mol.name+"-water")
    #    plt.plot(temperatures,h_i-h_f[(guest,structure)], label=structure+" "+mol.name+"-ice")
    plt.title("Latent heat per water molecule")
    plt.xlabel("Temperature [K]")
    plt.ylabel("\Delta H  [kJ/mol/water mol]")
    plt.legend()

    #per guest
    plt.figure(4)
    for guest,structure in targets:
        mol = moldict[guest]
        plt.plot(temperatures,dh_fwpg[(guest,structure)], label=structure+" "+mol.name+"-water")
        plt.text(temperatures[0],dh_fwpg[(guest,structure)][0], structure+" "+mol.name+"-water")
    plt.title("Latent heat per guest molecule")
    plt.xlabel("Temperature [K]")
    plt.ylabel("\Delta H  [kJ/mol/guest mol]")
    plt.legend()

    #per kg
    plt.figure(5)
    plt.plot(temperatures,(h_w-h_i)/18.0*1000, label="ice-water")
    plt.text(temperatures[0],(h_w[0]-h_i[0])/18.0*1000, "ice-water")
    for guest,structure in targets:
        mol = moldict[guest]
        plt.plot(temperatures,dh_fwpm[(guest,structure)], label=structure+" "+mol.name+"-water")
        plt.text(temperatures[0],dh_fwpm[(guest,structure)][0], structure+" "+mol.name+"-water")
    plt.title("Latent heat per mass")
    plt.xlabel("Temperature [K]")
    plt.ylabel("\Delta H  [kJ/mol/kG]")
    plt.legend()

    plt.show()


def LatentHeatMode():
    import interpolate
    #only one guest
    for guest in guests:
        #determine melting point and crystal type
        mp = 0
        mp_struc = ""
        for structure in structures:
            mp_w = interpolate.zerocross(temperatures, mu_f[(guest, structure)] - mu_w3)
            if debug:
              print mp_w,"MPW"
              print "mu_f",mu_f[(guest,structure)]
              print "mu_w3",mu_w3
              print "diff", mu_f[(guest,structure)] - mu_w3
            if temperatures[0] < mp_w < temperatures[-1]:
                if mp < mp_w:
                    mp_struc = structure
                    mp = mp_w
            # if temperatures[0] < mp_w < temperatures[-1]:
            #     if mp < mp_w and 273.15 < mp_w:
            #         mp_struc = structure
            #         mp = mp_w
            # mp_i = interpolate.zerocross(temperatures, mu_f[(guest, structure)] - mu_i)
            # if temperatures[0] < mp_i < temperatures[-1]:
            #     if mp < mp_i < 273.15:
            #         mp_struc = structure
            #         mp = mp_i
        if mp_struc != "":
          latent = interpolate.zerocross( dh_fwpm[(guest, mp_struc)], temperatures - mp )
          print("{0}\t{1}\t{2}\t{3}\t# kJ/kG".format(guest,mp_struc, mp,latent))



if mode == "Test":
  TestReport(guests)
  sys.exit(0)


if mode == "Graphs":
  GraphMode()

if mode == "Latent":
  LatentHeatMode()
