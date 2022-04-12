#coding: utf-8

import math
def relocate( box, water, tag, origin ):
    print "@BOX3"
    print box[0],box[1],box[2]
    print tag
    print len(water)
    for w in water:
        for dim in range(3):
            w[dim] -= origin[dim]
            w[dim] -= math.floor( w[dim] / box[dim] + 0.5 ) * box[dim]
        for x in w:
            print x,
        print

def usage():
    print "usage: %s waterconfig.nx4a < guestconfig.ar3a" % sys.argv[0]
    sys.exit(0)

import sys
if len(sys.argv) != 2:
    usage()
nx4a = open(sys.argv[1])
box = [0.0]*3
water = []
while True:
    line = nx4a.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if columns[0] == "@BOX3":
        line = nx4a.readline()
        box = map(float,line.split())
    elif columns[0] in ( "@NX4A", "@AR3A" ) :
        tag = columns[0]
        line = nx4a.readline()
        nmol = int(line)
        for i in range(nmol):
            line = nx4a.readline()
            columns = map(float,line.split())
            water += [columns]

while True:
    line = sys.stdin.readline()
    if len(line) == 0:
        break
    columns = line.split()
    if columns[0] == "@BOX3":
        line = sys.stdin.readline()
        box = map(float,line.split())
    elif columns[0] in ( "@NX4A", "@AR3A" ) :
        line = sys.stdin.readline()
        nmol = int(line)
        for i in range(nmol):
            line = sys.stdin.readline()
            columns = map(float,line.split())
            relocate( box, water, tag, columns[0:3] )
