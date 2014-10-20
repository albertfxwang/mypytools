#!/usr/bin/env python
#
""" Some image utilities """

__author__ = "Henry C. Ferguson, Space Telescope Science Institute"
__copyright__ = "Copyright 2012, Henry C. Ferguson"
__credits__ = ["Henry C. Ferguson"]
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Henry C. Ferguson"
__email__ = "ferguson@stsci.edu"
__status__ = "Development"

from numpy import *
from iterstat import *

def blocksum(m,k):
    """ Block sum over kxk pixels """
    return add.reduceat(add.reduceat(m,arange(0,m.shape[0],k),axis=0),arange(0,m.shape[1],k),axis=1)

def makemask(img,threshold=0.05):
    """ Mask all pixels above a given threshold """
    return img > threshold

def smoothsample_masked(img,mask,width=9, reject_sigma=3.0, niter=3):
    """ Put down non-overlapping square apertures over the entire unmasked portion of an image.
        Return a (possibly clipped) mean of the values in each aperture.
          - img: The image to measure (2-d array)
          - mask: The mask  (2-d array)
          - width: The width of the square apertures
          - reject_sigma: The threshold for iterative rejection
    """
    nx,ny = shape(img)
    w = width/2
    zz = []
    xx = []
    yy = []
    # Mask out the borders
    xg,yg = mgrid[0:nx,0:ny]
    mask = mask | (xg<w) | (xg>nx-w) | (yg<w) | (yg>ny-w)
    # Make a list of aperture centers for apertures that don't included masked pixels
    xg,yg = mgrid[0:nx/width,0:ny/width]
    mm = blocksum(mask,width)
    good = where(mm == 0)
    xlist = (xg[good]*width+width/2).tolist()
    ylist = (yg[good]*width+width/2).tolist()
    # Cycle through these apertures, measuring the mean flux
    if reject_sigma:
        for x,y in zip(xlist,ylist):
            zz += [iterstat(img[int(x-w):int(x+w)+1,int(y-w):int(y+w)+1],niter, 
                 low=reject_sigma, high=reject_sigma,justmean=True)]
    else:
        for x,y in zip(xlist,ylist):
            zz += [img[int(x-w):int(x+w),int(y-w):int(y+w)].mean()]
    xx = array(xlist)
    yy = array(ylist)
    zz = array(zz)
    return xx,yy,zz

