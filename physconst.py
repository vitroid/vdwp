
class physconst():
    def __init__(self):
        #constants
        #self.J2K = 0.120273 # K/J 20110521 modified.
        self.MassMe = 16e-3 # kg/mol

        #Physical propertires
        self.kB   = 1.380662e-23 # J/K
        self.NA   = 6.022045e23  # 1
        self.J2K =  1.0/(self.kB*self.NA) # K/J
        self.NkB  = self.kB * self.NA * 1e-3 # kJ/mol/K
        self.h    = 6.626176e-34 # J s Planck's h
        self.cc   = 2.99792e10


pc = physconst()

if __name__ == "__main__":
    from math import *
    T = 273.15
    print("a(Xe)",   -3 * pc.kB * T / 2 * log( (131.3*0.001/ pc.NA) * pc.kB * T * 2.0*pi/pc.h**2)*pc.NA*0.001) # kJ/mol
    print("a(Et)",   -3 * pc.kB * T / 2 * log( (30.0*0.001 / pc.NA) * pc.kB * T * 2.0*pi/pc.h**2)*pc.NA*0.001) # kJ/mol
    print("a(Pr)",   -3 * pc.kB * T / 2 * log( (44.0*0.001 / pc.NA) * pc.kB * T * 2.0*pi/pc.h**2)*pc.NA*0.001) # kJ/mol
    print("a(IsoBu)",-3 * pc.kB * T / 2 * log( (58.1*0.001 / pc.NA) * pc.kB * T * 2.0*pi/pc.h**2)*pc.NA*0.001) # kJ/mol
    print("a(CF4)",  -3 * pc.kB * T / 2 * log( (88.0*0.001 / pc.NA) * pc.kB * T * 2.0*pi/pc.h**2)*pc.NA*0.001) # kJ/mol
    print("c(Et)",(pc.kB*T*log(2)-pc.kB*T*log((30*0.001/pc.NA)*(1.530e-10)**2/4.0*pc.kB*T*2.0*pi/pc.h**2))*pc.NA*0.001) # kJ/mol
    #print (30*0.001/pc.NA)*(1.530e-10)**2/4.0
