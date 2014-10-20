#!/usr/bin/env python
import sys
from numpy import *
from gauss import *

def lognormal(x,mean,sigma):
    l = log10(x)
    mean = log10(mean)
    val = gauss(l,mean,sigma)
    return val

def random_lognormal(n,mean,sigma):
    mean = log10(mean)
    val = random.normal(l,mean,sigma)
    data = random.normal(size=n)
    data = data*sigma+mean
    data = exp(data)
    return data

if __name__ == "__main__":
    mean = float(sys.argv[1])
    sigma = float(sys.argv[2])
    r = arange(0.01,10.,0.01)
    y = lognormal(r,mean,sigma)
    for i in range(len(r)):
        print r[i],y[i]
    
