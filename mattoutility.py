from re import *
import math
import string
#from scipy import constants

spaces=compile(" +")

def spacesplit(line):
    line = string.rstrip(line," \t\r\n")
    columns=spaces.split(line)
    while columns[0] == "":
        columns.pop(0)
    return columns
