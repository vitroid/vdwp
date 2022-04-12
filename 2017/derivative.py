#!/usr/bin/env python


import numpy
from logging import getLogger

def derivatives(x,y):
    logger = getLogger()
    if len(x) != len(y):
        logger.warning("List size is different in " + self.__name__)
    n = len(y)
    #intermediate values
    d = numpy.zeros(n)
    d[1:n-1] = (y[2:n] - y[0:n-2]) / (x[2:n] - x[0:n-2])
    d[0]  = (y[1] - y[0]) / (x[1] - x[0])
    d[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
    return d
