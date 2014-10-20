#!/usr/bin/env python
# 
# Version 1.0.0 mlfit utilities that do not use scipy
# 4/01/09
# version 1.3.0
# Modified to use numpy by Kuang-Han Huang

""" 
Utility functions for maximum-likelihood fitting.
Updated from mlutil.py under $biv/galfit_transfer
"""

__author__ =  'Henry C. Ferguson / Kuang-Han Huang'
__version__=  '1.4.0'

from numpy import *
import numpy
import sys
import pyfits
import cPickle
from pygoods import *
import matplotlib as mpl
import matplotlib.pyplot as plt

ag = 1
try:
    from mpl_toolkits.axes_grid import AxesGrid
except:
    ag = 0
    pass

def writefits2d(a,fitsfile):
    """ writefits2d -- write out an array to a fits file
        a -- array (numarray or Numeric)
        fitsfile -- output filename
    """
    fitsobj = pyfits.HDUList()
    hdu = pyfits.PrimaryHDU()
    hdr = hdu.header
    hdr.update('crval1',0.)
    hdr.update('crval2',0.)
    hdr.update('crpix1',a.shape[0]/2.)
    hdr.update('crpix2',a.shape[1]/2.)
    
    # Check if we already have a numarray array, if not convert
    x = N.arange(2)
    print type(a)
    if type(a) != type(x):
        hdu.data = N.array(a.tolist())
    else:
        hdu.data = a.copy()
    fitsobj.append(hdu)
    fitsobj.writeto(fitsfile)
  

def gauss(x,m,sigma):
    sigma2 = sigma*sigma
    return exp(-(x-m)**2/(2.*sigma2))/sqrt(2.*pi*sigma2)

def limit_data1d(data,limits):
    flag = (data>min(limits)) & (data<max(limits))
    return N.compress(flag,data)

def limit_data2d(data,limits):
   """
   Prune the data points (given in a 2xN array) outside the limits.
   Returns another 2xN array.
   """
   # if 2D data: shape(data)=(2,ndata)
   if len(N.shape(data)) == 1:
      return limit_data1d(data,limits)
   crit = (data[0]>=limits[0,0])*(data[0]<limits[0,1])
   crit = crit * (data[1]>=limits[1,0])*(data[1]<limits[1,1])
   d0 = compress(crit,data[0])
   d1 = compress(crit,data[1])
   cdata = array([d0,d1])
   return cdata

def grid_data1d(inputdata,limits,nbins,verbose=0):
    """ Grid the data on a pixel grid.

        Arguments: 
        limits - a list of min,max 
        nbins - the number of bins
    """
    data = limit_data1d(inputdata,limits)
    ndata = len(data)
    coords = N.zeros(len(data))
    output = N.zeros(nbins)
    xmin = limits[0]
    xmax = limits[1]
    dx = (xmax-xmin)/nbins
    coords = ((data-xmin)/dx).astype('int')
    for i in range(len(data)):
       output[coords[i]] += 1 
    if verbose:
        print "ndata,len(output) ",ndata,len(output)
    return (output,ndata)

def grid_data2d(inputdata,limits,pixdx):
    """ 
    Grid the data on the same pixel grid as the models.
    inputdata -- a 2xN array
    limits -- the boundaries of the parameter space (including the upper limits ine each dimension)
    pixdx -- the "pixel" size in each dimension, a length-2 array
    """
    data = limit_data2d(inputdata,limits)
    xpix = arange(limits[0,0],limits[0,1]+pixdx[0],pixdx[0])  # edges of pixels in the 1st dimension, including upper edge of the last pixel
    ypix = arange(limits[1,0],limits[1,1]+pixdx[1],pixdx[1])
    datagrid = histogram2d(inputdata[0],inputdata[1],bins=[xpix,ypix])[0]
    return datagrid
        
def normalize(model,value):
    """ Return a model where the integral in the unmasked region is
        equal to the total number of data points. 

        Arguments:
        ndata - number of data points in the masked region
        model - masked array of unnormalized model values 
    """
    total = sum(model.ravel())
    return model*value/total
        
def poissonprob(n,mu):
    """ Return Poisson probability of n events given an expectation of mu.
        For mu > 50 return Gaussian probability.
    """
    maxP = 100.
    gmask = n >= maxP
    npoisson = choose(gmask,(n,maxP))
    # Poisson estimates
    poissonresult = (mu**n * exp(-mu))/factorial(npoisson)
    # Gaussian estimates
    ngauss = choose(gmask,(maxP,n))
    gaussresult = gauss(mu,ngauss,sqrt(mu))
    result = choose(gmask,(poissonresult,gaussresult))
    return result

def poissonlr(data,model):  # keep?
    """Return the -sum(log(poisson probability))
       Dolphin, 2002, MNRAS 332, 91, equation 10.

       Arguments:
       data - masked data array (integers -- number of points per bin)
       inputmodel - masked model array (pdf normalized to same total number) 
       Note: this will below up if any bin in model == 0
    """
    i = nonzero(data.ravel())
    d = take(data.ravel(),i)
    m = take(model.ravel(),i) 
    log_plr = m-d+d*log(d/m)
    total_log_plr = sum(2.*log_plr.ravel())
    return total_log_plr  

def poissonlikelihood(data,inputmodel):
    """Return the -sum(log(poisson probability))

       Arguments:
       data - masked data array (integers -- number of points per bin)
       inputmodel - masked model array (pdf normalized to same total number) 
       Note: this will below up if any bin in model == 0
    """
    l = log(poissonprob(data,inputmodel))
    total_logl = sum(l.ravel())
    return -total_logl

def draw_from_pdf(ntot,model,limits,pixdx,verbose=0):
   """ 
   Return ntot data points drawn from a model distribution function. This is general for n-d model arrays.

   Arguments:
   ntot -- number of data values to create
   model -- a N-d array of PDF from which to draw points
   limits -- list of upper & lower limits for data, in the format of [[xlo,xhi],[ylo,yhi],[zlo,zhi],...]
   pixdx -- the pixel size of the model in each dimension
   Note that model, limits, and pixdx need to be consistent with each other, i.e.
   shape(model) == (limits[:,1]-limits[:,0])/pixdx
   """
   xmin = limits[:,0]
   xmax = limits[:,1]
   pdf = model.copy()
   pdf = pdf/max(pdf.ravel())    # set largest value to 1.0
   dimensions = shape(pdf)
   print "dimensions", dimensions
   ndimensions = len(dimensions)	# Number of dimensions
   datapoints = N.zeros((ndimensions,ntot),"float")
   index = N.zeros(ndimensions,"int")
   n = 0
   for n in range(ntot):                # Populate with random draws
      drawn = 0
      while drawn == 0:
         # randomly draw a point within the model boundaries, from a flat distribution (could also implement other distributions later?)
         x = numpy.random.uniform(limits[:,0],limits[:,1])
         index = floor((x-limits[:,0])/pixdx).astype('int')
         index_ravel = N.ravel_multi_index(index,shape(pdf))  # convert to the index in the raveled 1-D array 
         y = numpy.random.uniform()   # randomly draw between 0 and 1
         if pdf.ravel()[index_ravel] > y:  # then draw a point from this bin
            datapoints[:,n] = x
            drawn = 1
   return datapoints

def findindex2d(limits,pixdx,x0):
    """
    Find the corresponding index of the pixels containing x0,
    and the coordinates (lower edges) of these pixels.
    x0 should be a Nx2 array, where N is the number of points being queried.
    returns the indices and coordinates of the pixels (all are Nx2 arrays)
    Note that if any x0 is outside of limits, the returned indices will be coerced to the 
    those of the pixels at the boundary.
    """
    x0 = array(x0)
    id = (x0 - limits[:,0])/pixdx
    id = floor(id).astype('int'); print id
    xcoord = arange(limits[0,0],limits[0,1],pixdx[0])  # the x-coordinates of all pixels
    ycoord = arange(limits[1,0],limits[1,1],pixdx[1])  # the y-coordinates of all pixels
    #coord = array([limits[0,0]+pixdx[0]*id[0],limits[1,0]+pixdx[1]*id[1]])
    coord = array([xcoord[id[:,0]],ycoord[id[:,1]]]).swapaxes(0,1)
    # coerce all pixel indices outside of the limits to the end points in each dimension
    coord = minimum(coord,[xcoord[-1],ycoord[-1]])
    coord = maximum(coord,[xcoord[0],ycoord[0]])
    return id, coord

def testfunc((x,y,z),power,offset):
    return (x+y*z)**power+offset

def bincoords2d(limits,cell_widths):
   """ 
   Create arrays of cell coordinates.
   Arguments:
   limits -- array of limits (min,max); contains the upper limit of the last bin
   cell_widths -- the width of each cell in both dimensions
   Returns:
      x0 -- lower param-space coordinates of each cell
      x1 -- upper param-space coordinates of each cell
      xi -- integer indices of each cell
   x0, x1, xi are all (nobj, 2) numpy arrays.
   """
   # Create an array of bin coordinates
   # xi are indices, x0 and x1 are the lower and upper boundary of each bin
   XX = arange(limits[0,0],limits[0,1],cell_widths[0])
   YY = arange(limits[1,0],limits[1,1],cell_widths[1])
   ZZ = meshgrid(XX,YY)
   xxi = arange(len(XX))
   yyi = arange(len(YY))
   x0 = column_stack((ZZ[0].ravel(),ZZ[1].ravel()))
   x1 = x0 + cell_widths
   xi = meshgrid(xxi,yyi)
   xi = column_stack((xi[0].ravel(),xi[1].ravel()))
   
   return (x0,x1,xi)

def rlognormal(rmin,rmax,mag,(sigma,r0,mag0,beta)):
    lumratio = 10.**((mag-mag0)/-2.5)
    mu = log(r0)+beta*log(lumratio)
    val = random.normal(mu,sigma)
    radius = exp(val)
    return radius

def roundid(lowlim,uplim,step,value):
    g = N.arange(lowlim,uplim+step,step)
    x = (value-lowlim)/step
    near = int(round(x))
    if value >= g[near]:
        id = near
    else:
        id = near - 1
    return id

if ag:
    def display_kernels(kgrid,n1,n2):
        # default: n1=12(mag), n2=17(Re)
        fig = figure(1,figsize=(12,8))
        grid = AxesGrid(fig,(0.1,0.15,0.9,0.58),
                        nrows_ncols = (n2,n1),
                        axes_pad = 0.0,
                        share_all=True,aspect=False)
                        
        #grid.set_aspect(0.6)
        for i in arange(n1):
            for j in range(n2)[::-1]:
                kernel = maximum(kgrid.kernels[(i,j)].kernel,1.e-50)
                if kgrid.kernels[(i,j)].completeness<0:
                    kernel = kernel * 0.
                else:
                    #nf = 1.0/max(kernel.ravel())
                    kernel = kernel * kgrid.kernels[(i,j)].completeness
                #grid[(n2-j-1)*n1+i].
                grid[(n2-j-1)*n1+i].imshow(kernel.swapaxes(0,1),origin='lower',
                                           vmax=0.002,aspect=0.4)
                #grid[(n2-j-1)*n1+i].xscale=1.5
                #grid[(n2-j-1)*n1+i].yscale=1.0
                grid[(n2-j-1)*n1+i].get_aspect()
        
        grid.axes_llc.set_xticks([])
        grid.axes_llc.set_yticks([])
        draw()
       
