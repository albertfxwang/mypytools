
# radial coordinate = (r-L)/(r+L)
# 

from numpy import *
from numpy.polynomial import chebyshev as c
from scipy.ndimage import interpolation
from scipy import stats
from pygoods import *

def TL(r,L,coeffs):
    print "coeffs = ",coeffs
    x = (r-L)/(r+L)
    return c.chebval(x,coeffs)

def chef(r,theta,L,fs,fc):
    result = 0
    for m in range(fs.shape[1]):
       print "m,fc,fs = ", m, fc[:,m], fs[:,m]     
       term1 = TL(r,L,fc[:,m])*cos(m*theta)
       term2 = TL(r,L,fs[:,m])*sin(m*theta)
       result = result+term1+term2
    return result

def tstchef(fs,fc):
    L = 3.
    xc = 50.
    yc = 50.
    x,y = mgrid[0:100,0:100]
    r = sqrt((x-xc)**2+(y-yc)**2)
    theta = arctan2(y-yc,x-xc)
    return chef(r,theta,L,fs,fc)

def testimage1(N=500):
    x,y = mgrid[0:N,0:N]
    gauss = stats.norm(250.,scale=4.).pdf(x)*stats.norm(250.,scale=8.).pdf(y)
#   out = zeros(gauss.shape)
#   interpolation.rotate(gauss,45.,output=out,reshape=False)
    out = gauss
    return out
