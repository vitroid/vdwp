#!/usr/bin/env python


def loadHIST(file):
    line = file.readline()
    dimen = int(line)
    line = file.readline()
    columns = line.split()
    nbin = int(columns[0])
    vmin, binwidth = [float(x) for x in columns[1:3]]
    histo = dict()
    for i in range(nbin):
        line = file.readline()
        weight = float(line)
        value = vmin + (i+0.5)*binwidth
        histo[value] = weight
    return histo


def loadSHST(file):
    line = file.readline()
    dimen = int(line)
    line = file.readline()
    columns = line.split()
    vmin, binwidth = [float(x) for x in columns[0:2]]
    histo = dict()
    while True:
        line = file.readline()
        columns = line.split()
        i = int(columns[0])
        if i < 0:
            break
        weight = float(columns[1])
        value = vmin + (i+0.5)*binwidth
        histo[value] = weight
    return histo


def loadAHisto(file):
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        columns = line.split()
        if len(columns) > 0:
            if columns[0] == "@HIST":
                return loadHIST(file)
            elif columns[0] == "@SHST":
                return loadSHST(file)
    return None



def test():
    h1 = loadAHisto(open("data/LJPR2LB_.16hedra.histo"))
    for key in sorted(h1.keys()):
        print(key,h1[key])



if __name__ == "__main__":
    test()
