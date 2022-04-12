#!/usr/bin/env python

import numpy


#assume the function increases monotonically
def zerocross(X,Y):
    if Y[0] > 0:
        return X[0] - 1
    if Y[-1] < 0:
        return X[-1] + 1
    for i in range(0,len(Y)-1):
        z = Y[i]*Y[i+1]
        if z <= 0:
            delta = X[i+1]-X[i]
            slope = (Y[i+1]-Y[i]) / delta
            dx = -Y[i] / slope
            return X[i] + dx
    print("ERROR zerocross",X,Y)


#find the value at x
def value(x,X,Y):
    return zerocross(Y,X-x)


#X = numpy.array([1,2,3,4,5])
#Y = numpy.array([-3,-1,-4,-1,-5])
#print(value(2.5,X,Y))

