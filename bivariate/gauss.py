#!/usr/bin/env python
from numpy import *

def gauss(x,m,sigma):
    sigma2 = sigma*sigma
    return exp(-(x-m)**2/(2.*sigma2))/sqrt(2.*pi*sigma2)

def normgauss(x,m,sigma):
    g = gauss(x,m,sigma)
    return g/sum(g)
