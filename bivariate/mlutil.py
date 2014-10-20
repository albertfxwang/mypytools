#!/usr/bin/env python
# 
# Version 1.0.0 mlfit utilities that do not use scipy
# 4/01/09
# version 1.3.0
# Modified to use numpy by Kuang-Han Huang

""" Utility functions for maximum-likelihood fitting.
"""

__author__ =  'Henry C. Ferguson'
__version__=  '1.3.0'

from numpy import *
#import MA 
#import RandomArray
#import random
import numpy as N
import sys
import pyfits
import cPickle
from pygoods import *
from pylab import *

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
  
    
#def tonumarray(a):
#    return array(buffer(a),type=a.dtype.char,shape=a.shape)

#def tonumarrays(*arraylist):
#    ret = []
#    for a in arraylist:
#        ret= ret + [array(buffer(a),type=a.typecode(),shape=a.shape)]
#    return ret

def gauss(x,m,sigma):
    sigma2 = sigma*sigma
    return exp(-(x-m)**2/(2.*sigma2))/sqrt(2.*pi*sigma2)

def limit_data1d(data,limits):
    flag = (data>min(limits)) & (data<max(limits))
    return N.compress(flag,data)

def limit_data(data,limits):
    # if 2D data: shape(data)=(2,ndata)
    if len(N.shape(data)) == 1:
        return limit_data1d(data,limits)
    crit = (data[0]>=limits[0,0])*(data[0]<limits[0,1])
    crit = crit * (data[1]>=limits[1,0])*(data[1]<limits[1,1])
    d0 = compress(crit,data[0])
    d1 = compress(crit,data[1])
    cdata = array([d0,d1])
    return cdata
    #cdata=[]
    #flag = N.array([True]).repeat(shape(data)[1])
    #for i in range(len(data)):
    #    xmin = min(limits[i])
    #    xmax = max(limits[i])
    #    iflag = (data[i]>xmin) & (data[i]<xmax)
    #    flag = flag*iflag
    #ndata = sum(flag)
    #for i in range(len(data)):
    #    cdata += [N.compress(flag,data[i])]
    #return N.array(cdata).reshape(len(data),ndata)
    
def mask_data(data,boundaries):
    """ Return a data array that has been purged of values outside the boundaries. 
 
        Arguments:
        boundaries(data) -- a user-supplied function that takes
        the data array and returns a mask (1 = accept, 0=reject)
    """
    b = boundaries(data)
    return N.compress(b,data,dimension=1) 

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

def grid_data(inputdata,limits,pixdx):
    """ Grid the data on the same pixel grid as the models.
    """
    data = limit_data(inputdata,limits)
    npix = (limits[:,1]-limits[:,0])/pixdx
    datagrid = N.zeros(npix,"int")
    data = data.swapaxes(0,1)   # reshape shape(data) to (ndata,2)
    for i in range(len(data)):
        ind,c = findindex(limits,pixdx,data[i])
        datagrid[ind[0],ind[1]] += 1
    return datagrid
    

def grid_data_old(inputdata,limits,nbins):
    """ Grid the data on the same pixel grid as the models.
 
        Arguments:
        grid - a list of min,max,nbin for each coordinate
        nbins - a list of the number of bins in each axis
    """
    # Convert limits, nbins to arrays if data is one-dimensional
    if len(N.shape(inputdata)) == 1:
        return (grid_data1d(inputdata,limits,nbins))
    data = limit_data(inputdata,limits)
    ndata = N.shape(data)[1]
    coords = N.zeros(N.shape(data),dtype='int')
    output = {}
    for i in range(len(data)):
        xmin = limits[i,0]
        xmax = limits[i,1]
        nbin = nbins[i]
        dx = (xmax-xmin)/nbin
        coords[i] = ((data[i]-xmin)/dx).astype('int')
    # Now probably there is a clever way to do the next bit without
    # a loop, but I can't think of it. For now, we just step through
    # the data points and increment the appropriate value in the grid
    # by one for each data point.
    coords = coords.swapaxes(0,1)
    for i in range(N.shape(data)[1]):
        kc = tuple(coords[i])
        if not output.has_key(kc):
            output[kc] = 0
        output[kc] += 1 
    return (output,ndata)
        
def normalize(model,ndata):
    """ Return a model where the integral in the unmasked region is
        equal to the total number of data points. 

        Arguments:
        ndata - number of data points in the masked region
        model - masked array of unnormalized model values 
    """
    total = sum(model.ravel())
    return model*ndata/total
        
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

def poissonlr(data,model):
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

def draw_from_pdf(ntot,model,limits,verbose=0):
    """ Return ntot data points drawn from a model distribution function.

        Arguments:
        ntot -- number of data values to create
        model -- a bivmodel class object
        limits -- list of upper & lower limits for data
        nbins -- number of bins
        model must already be gridded on the grid given by limits,nbins
    """
    limits = N.array(limits)
    modlimits = model.limits
    dx = model.pixdx
    xmin = limits[:,0]
    xmax = limits[:,1]
    idmin = (xmin-modlimits[:,0])/dx
    idmin = idmin.astype('int')
    idmax = (xmax-modlimits[:,0])/dx
    idmax = idmax.astype('int')
    npix = shape(model.model)   # model is super-imposed on the parameter space
    #dx = (xmax-xmin)/npix              # calculate the grid spacing
    #pdf = model.model[idmin[0]:idmax[0]+1,idmin[1]:idmax[1]+1]
    pdf = model.model[idmin[0]:idmax[0],idmin[1]:idmax[1]]
    pdf = pdf/max(pdf.ravel())    # set largest value to 1.0
    dimensions = shape(pdf)
    print "dimensions", dimensions
    ndimensions = len(dimensions)	# Number of dimensions
    #pdf = model.model/max(model.ravel())         # Set largest value to 1
    # draw_from_pdf only depends on the relative values in the model, is
    # indifferent to overall normalization
    datapoints = N.zeros((ndimensions,ntot),"float")
    index = N.zeros(ndimensions,"int")
    n = 0
    for n in range(ntot):                # Populate with random draws
        drawn = 0
        while drawn == 0:
            for i in range(ndimensions):
                # randomly choose a pixel for consideration
                index[i] = int(round(N.random.uniform(0,dimensions[i]-1)))   
                # the index of the pixel we're looking at
            y = N.random.random()  
            if pdf[index[0],index[1]] > y:  # then draw a point from this bin
                epsilon = dx*N.random.random(size=len(index)) 
                # randomly generate an offset from pixel boundary
                d = index*dx + xmin + epsilon*dx  
                # computes the value of a random point in the bin
                datapoints[:,n] = d
                drawn = 1
    return datapoints

def findindex(limits,pixdx,x0):
    """Find the corresponding index of coord x0"""
    notwithin = (x0<limits[:,0])|(x0>limits[:,1])
    if any(notwithin):
        raise ValueError, "cannot find index"
    id = (x0 - limits[:,0])/pixdx
    id = floor(id).astype('int')
    coord = array([limits[0,0]+pixdx[0]*id[0],limits[1,0]+pixdx[1]*id[1]])
    return id, coord
    
def readkgrid(kgridfile):
    f = open(kgridfile)
    kgrid = cPickle.load(f)
    f.close()
    return kgrid

def testfunc((x,y,z),power,offset):
    return (x+y*z)**power+offset


class bivmodel:
    """An object that contains information of model limits, pixel size, etc.
    """
    def __init__(self,limits,pixdx,model_array):
        self.limits = limits
        self.pixdx = pixdx
        self.model = model_array
        mags = N.arange(limits[0,0],limits[0,1],pixdx[0])
        logr = N.arange(limits[1,0],limits[1,1],pixdx[1])
        coords = N.zeros((len(mags),len(logr),2))
        for i in range(len(mags)):
            coords[i,:,0] = mags[i]
        for j in range(len(logr)):
            coords[:,j,1] = logr[j]
        self.coords = coords  
        # self.coords[i,j] returns (mag,logr) corresponding to the (i,j)th
        # pixel
    def copy(self):
        x = self.__init__(self.limits,self.pixdx,self.model)
        return x

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

def plotcatalog(cat,limits,pixdx,color='black'):
	# plot the GALFIT catalog points
    s = sextractor(cat)
    mag = s.mag
    magauto = s.mag_auto
    re = s.re
    resex = s.flux_radius_1
    crit = (mag>=limits[0,0])*(mag<limits[0,1])*(log10(re)>=limits[1,0])*(log10(re)<limits[1,1])
    crit = crit*(s.chisqnu<=0.45) * (abs(s.mag_err)<=1.0)
    crit = crit * ((s.cullflag==0)+(s.cullflag==11)+(s.cullflag==14)) # Only the best ones after visual culling 
    crit = crit * (abs(s.re_err/s.re)<=1.0)
    #crit = crit * (s.sersicn>0.1) #* (s.sersicn<7.95)
    mag = N.compress(crit,mag)
    logre = N.compress(crit,log10(re))
    magauto = N.compress(crit,magauto)
    logresex = N.compress(crit,log10(resex))
    # plot the SE parameters of GALFIT catalog points
    print "len(mag)",len(mag)
    mag_i = (mag-limits[0,0])/pixdx[0]
    logre_i = (logre-limits[1,0])/pixdx[1]
    magauto_i = (magauto-limits[0,0])/pixdx[0]
    logresex_i = (logresex-limits[1,0])/pixdx[1]
    p1=plot(mag_i,logre_i,'.',color=color,markersize=2) #ffff33
    #p1=plot(magauto_i,logresex_i,'.',color='black',markersize=4)
    xlim(0,275)
    ylim(0,150)
    #xticks(arange(0,300,50),arange(21.0,27.0))
    #yticks([50,100,150],[0.03,0.3,3.0])
    return p1

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
       
