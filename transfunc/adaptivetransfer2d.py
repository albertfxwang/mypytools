#!/usr/bin/env python
# 
# From a list of model->observed values (e.g. results of artifical star 
# tests) construct convolution kernels. 
#
# Input values are assumed to be uniformly distributed.
#
# 2/08/05
# 
# Started from v0.8.0 of transfer.py, and added adaptive kernel
# for 2-d case.
#
# 9/7/10
# version 2.0
# Re-written by Kuang-Han Huang for numpy

from numpy import *
import KPDFadaptnumpy as KPDF
import sys
import scipy
import draw_ndist as dn
#from uniparam import limits,nbins,pixdx
from pygoods import *
import copy, cPickle

__author__ =  'Henry C. Ferguson/Kuang-Han Huang'
__version__=  '2.0' 

def normalize(model,ndata):
    """ Return a model where the integral in the unmasked region is
        equal to the total number of data points. 

        Arguments:
        ndata - number of data points in the masked region
        model - masked array of unnormalized model values 
    """
    total = sum(model.ravel())
    return model*ndata/total

def getdx(limits,nbins):
    """Could be used to return bin width or pixel size"""
    return (limits[:,1]-limits[:,0])/nbins

class kernelgrid:
    """ An evenly spaced grid of convolution kernels """
    def __init__(self,limits,nbins,pixdx):
        self.limits=limits
        self.nbins = nbins
        self.dxgrid = getdx(limits,nbins)
        self.pixdx = pixdx
        self.kernels = {}
    def getkernel(self,x):  
        """ Return the kernel corresponding to a given set of coords.
        """
        y = (x-self.limits[:,0]) / self.dxgrid
        i = y.astype("int")
        if self.kernels.has_key(tuple(i)):
            return self.kernels[tuple(i)]
        else:
            print "No kernel for ",x
            return -99

class kernel:
    """ A convolution kernel associated with a certain range of parameter space
    """
    def __init__(self,x0_bin,x1_bin,pixdx,completeness,kernel_array):
        self.x0_bin = x0_bin                         # Where to start using this kernel
        self.x1_bin = x1_bin                         # Where to stop using this kernel
        self.pixdx = pixdx                         # kernel pixel size 
        self.i0 = -array(shape(kernel_array))/2
        self.i1 = array(shape(kernel_array))/2
        self.rx0 = -array(shape(kernel_array))*pixdx/2.
        self.rx1 = array(shape(kernel_array))*pixdx/2.
        self.completeness = completeness
        self.kernel = kernel_array
    def update(self,kernel_array):
        s1 = array(shape(self.kernel)).astype("float")
        s2 = array(shape(kernel_array)).astype("float")
        self.pixdx = self.pixdx*s1/s2
        self.i0 = -array(shape(kernel_array))/2
        self.i1 = array(shape(kernel_array))/2
        self.rx0 = -array(shape(kernel_array))*self.pixdx/2.
        self.rx1 = array(shape(kernel_array))*self.pixdx/2.
        self.kernel = kernel_array

def bincoords(limits,nbins):
    """ Create arrays of bin coordinates.
        Arguments:
        limits -- array of limits (min,max)
        nbins -- a 1-d numpy array of scalar number of bins in each dimension
    """
    # Create an array of bin coordinates
    # xi are indices, x0 and x1 are the lower and upper boundary of each bin
    dxgrid = getdx(limits,nbins)
    xi = indices(nbins.tolist())*1.0
    x0 = zeros(shape(xi),"float")
    x1 = zeros(shape(xi),"float")
    for i in range(len(xi)):
        x0[i] = xi[i]*dxgrid[i]+limits[i,0]
        x1[i] = xi[i]*dxgrid[i]+limits[i,0]+dxgrid[i]
    return (x0,x1,xi)

def iterstats(x,lowsig,highsig,niter):
    mean = scipy.stats.mean(x)
    sdev = scipy.stats.samplestd(x)
    for i in range(niter-1):
       select = ((x-mean) < highsig) & ((mean-x) < lowsig)
       x = compress(select,x) 
       mean = scipy.stats.mean(x)
       sdev = scipy.stats.samplestd(x)
    return (mean,sdev)

class transfer_grid_2d(Ftable):
    def __init__(self,simcat,limits,binwd,pixdx,xmax=array([3.0,1.0]),verbose=1,
        kernel_estimator="gaussian",chinu_max=0.4,reerrr_max=0.6,
        diskfrac=0.7,ntot=-1,sersicnbins=arange(0.,9.),method='galfit'):
        self.catname=simcat
        self.limits=limits  # the limits in each dimension
        self.binwd=binwd  # an array of bin widths
        self.pixdx=pixdx  # "pixel" size
        self.xmax=xmax    # the extent of each kernel
        self.verbose=verbose
        self.kernel_estimator=kernel_estimator
        self.chinu_max=chinu_max
        self.reerrr_max=reerrr_max
        self.diskfrac=diskfrac
        self.ntot=ntot  # ntot==-1 means use ALL points in the simulation catalog
        self.sersicnbins=sersicnbins
        self.method=method
        Ftable.__init__(self,simcat)  # read in the catalog

    def __deepcopy__(self):
        self2 = __init__(self,self.simcat,limits,binwd,pixdx,xmax=array([3.0,1.0]),verbose=1,
            kernel_estimator="gaussian",chinu_max=0.4,reerrr_max=0.6,
            diskfrac=0.7,ntot=-1,sersicnbins=arange(0.,9.),method='galfit')

    def __call__(self):
        """
    Function: takes simulation input & output and construct 2-D transfer function
    kernel grid

    **Updated for numpy**
    Arguments:
    simcat -- simulation output catalog
    limits -- limits within which to construct kernels -- 2-d array
    binwd -- width of each bin 
    pixdx -- pixel size of kernels 
    xmax -- maximum width of kernel
    chinu_max -- maximum chisq_nu for good qflags
    reerrrmax -- maximum (re_err/reout) for good qflags
    diskfrac -- disk fraction in simulation input
    ntot -- total number of simulation points drawn to make kernels
     
    Notes:
    completeness in each bin is (# of recovered & good for use objects)
       / (# of input objects in this bin). So this accounts for (1)
       SE detection completeness, (2)GALFIT recovery completeness, and
       (3)GALFIT measurement "quality completeness"
       Does NOT include any kind of sample selection (e.g. LBG 
       selection) completeness
        """
        nbins = around((self.limits[:,1]-self.limits[:,0])/self.binwd).astype('int')
        print "Kernel estimator form:", self.kernel_estimator
        print "Kernel grid limits:", self.limits
        print "Number of kernel:", nbins
        print "Kernel width:", self.xmax 

        last_invcovariance = ones(shape(self.limits),"float")

        # Draw points from simulation to match Sersic n distribution of LBGs
        # Then clean up using GALFIT quality flags
        gtype = self.galaxy_type
        nout = self.n_out
        if self.ntot<0:
            outindex = arange(len(self.d))
        else:
            outindex = dn.draw_ndist(self.diskfrac,self.ntot,gtype,nout,nbins=self.sersicnbins)

        # update arrays of drawn points
        mag_in = self.mag_in.take(outindex)
        mag_out = self.mag_out.take(outindex)
        mag_auto = self.mag_auto.take(outindex)
        #magerr = self.magouterr.take(outindex)
        re_in = self.re_in.take(outindex)
        re_out = self.re_out.take(outindex)  # GALFIT-measured effective radius
        re_se = self.flux_radius_1.take(outindex)  # SExtractor-measured half-light radius
        re_out_err = self.re_out_err.take(outindex)
        recovered_gf = self.recovered_galfit.take(outindex)
        recovered_se = self.recovered_se.take(outindex)
        if hasattr(self,'chisqnu'):
            chisqnu = self.chisqnu.take(outindex)
        npoints = len(mag_in)  # should equal to ntot
        siminput = array([mag_in,log10(re_in)])
        # GALFIT output
        print self.method
        if self.method == 'galfit':
            simoutput = array([mag_out,log10(re_out)])   
            residuals = array([mag_out-mag_in,log10(re_out)-log10(re_in)])
        # SExtractor output
        if self.method == 'sextractor':
            simoutput = array([mag_auto,log10(re_se)])
            residuals = array([mag_auto-magin,log10(re_se)-log10(re_in)])

        # Make GALFIT quality flags
        qflags = ((re_out_err/re_out)<=self.reerrr_max)
        if hasattr(self,'chisqnu'):
            qflags = qflags & (chisqnu<=self.chinu_max)
        qflags = qflags & (re_out>0) & (re_out_err!=inf)
        if hasattr(self,'x_flag'):  # if using a newer version of GALFIT
            qflags = qflags & (self.re_flag==0) & (self.mag_flag==0) & (self.n_flag==0)

        ## INITIALIZE KERNEL GRID

        dxgrid = getdx(self.limits,nbins)  #dxgrid - bin width in each dimension
        kgrid = kernelgrid(self.limits,nbins,self.pixdx)
        kgrid.kernel_estimator = self.kernel_estimator
        kgrid.chinu_max = self.chinu_max
        kgrid.reerrmax = self.reerrr_max
        kernel.diskfrac = self.diskfrac
        kmaxbins = around(self.xmax/self.pixdx).astype('int')
    
        # Create an array of bin coordinates in the kernel grid
        x0,x1,xi = bincoords(self.limits,nbins)   # x0,x1,xi will be 3-d arrays
        if self.verbose > 1:
            print "x0 = ",x0
            print "x1 = ",x1
            print "xi = ",xi
 
        # Create list of bin ranges and indices in the array of kernels
        # fx0,fx1,fxi will be 2-d arrays
        # fx0[:,i] will be the coordinates of the lower-left point of the ith grid
        # fx1[:,i] will be the coord of the upper-right point of the ith grid
        # fxi[:,i] will be the index of the ith grid
        fx0 = reshape(x0,(shape(x0)[0],len(x0.ravel())/shape(x0)[0]))
        fx1 = reshape(x1,(shape(x1)[0],len(x1.ravel())/shape(x1)[0]))
        fxi = reshape(xi,(shape(xi)[0],len(xi.ravel())/shape(xi)[0]))

        # Iterate through each bin, calculate completeness, 
        # clean up using GALFIT quality flags, then make kernels

        for i in range(shape(fx0)[-1]):
            # setup bin boundaries and indices
            bx0_bin = fx0[:,i].copy()   # Lower bound for data used to make this kernel
            bx1_bin = fx1[:,i].copy()   # Upper bound
            bx0_sw = fx0[:,i].copy()  # Copies if bin boundaries are expanded
            bx1_sw = fx1[:,i].copy()
            bxi = fxi[:,i].astype('int')

            binx0_bin = outer(bx0_bin,ones(npoints))
            binx1_bin = outer(bx1_bin,ones(npoints))

            # determines which input points are in the given bin
            withinbin = (siminput >= binx0_bin) & (siminput < binx1_bin)  # get data within bin
            withinbin = multiply.reduce(withinbin)  # multiply every element in withinbin together
            print "sum(withinbin)", sum(withinbin)

            # CALCULATE COMPLETENESS
            # completeness = (num. of objects detected by SE & pass GALFIT quality check) / 
            # (total num. of objects in that bin)
            n_recovered_gf_bin = sum(withinbin & recovered_gf)
            n_recovered_se_bin = sum(withinbin & recovered_se)
            n_goodfit_gf_bin = sum(withinbin * recovered_gf * qflags)
            print "Number recovered by SExtractor in this bin:", n_recovered_se_bin
            print "Number recovered by GALFIT in this bin:", n_recovered_gf_bin
            print "Number pass GALFIT quality check in this bin:", n_goodfit_gf_bin
            if n_goodfit_gf_bin > 4:
                #completeness = float(n_goodfit_gf_bin) / float(sum(withinbin)) 
                # including SE detection incompleteness
                completeness = float(n_goodfit_gf_bin) / float(n_recovered_se_bin) 
                # NOT counting SE detection incompleteness
            else:
                completeness = 0.0
            print "completeness,",completeness

            # retain points within this bin that are recovered by GALFIT
            # and pass quality flag test
            binflags = withinbin & recovered_gf & qflags
            bin_npoints = sum(withinbin & recovered_gf) 
            # within bin && recovered_gf && well-fit by GALFIT
            m0_bin = compress(binflags,siminput[0,:])                   
            m1_bin = compress(binflags,siminput[1,:])
            r0_bin = compress(binflags,residuals[0,:])
            r1_bin = compress(binflags,residuals[1,:])
            m_bin = array([m0_bin,m1_bin])
            r_bin = array([r0_bin,r1_bin])
            #rcflag_bin = compress(withinbin,recoveredflag) # 1=recovered, 0=not recovered
            nmodel_bin = len(r0_bin)   # number of points in the bin used to make kernel 

            # INITIALIZE KERNEL PROPERTIES
            karray = zeros(kmaxbins,"float")           # Create an empty kernel
            center = self.xmax/2. # the index of center of kernel
            # e.g. if kernel is 128 by 128 array, then index of center will be [64,64]
            kx = arange(-self.xmax[0]/2.+self.pixdx[0]/2.,self.xmax[0]/2.+self.pixdx[0]/2.,self.pixdx[0])
            ky = arange(-self.xmax[1]/2.+self.pixdx[1]/2.,self.xmax[1]/2.+self.pixdx[1]/2.,self.pixdx[1])
            gdata = KPDF.MPDF2DGrid2Array(kx,ky) # output PDF grid
            print "shape(gdata)",shape(gdata)
            covariance=cov(r_bin,rowvar=1)

            if n_goodfit_gf_bin <= 4:  # less than 4 points in this bin
                print "*** WARNING -- kernel has less than 4 points***"
                print "bx0_bin, bx1_bin = ",bx0_bin,bx1_bin
                pdf2d = ones((shape(kx)[0],shape(ky)[0]),"float")
                print "shape(pdf2d)", shape(pdf2d)
                lam = None
                covariance = None
            else: 
                rr = transpose(r_bin)
                if self.verbose > 1:
                    print "kx=",kx
                #covariance = cov(rr)

                # Compute scale factors to use for KPDF 
                # This is a kludge to prevent crashing on singular covariance matrices
                try:  
                    invcovariance = linalg.inv(covariance)
                except: 
                    print "*** Warning, singular inverse convariance matrix, trying Moore-Penrose inverse ***" 
                    try:
                        invcovariance = linalg.pinv(covariance,
                                rcond=1.e-10)
                    except:
                        # This will bomb if it happens on the first pass
                        print "*** No joy...using last good invcovariance ***" 
                        invcovariance = last_invcovariance
                last_invcovariance = invcovariance
          
                renorm_constant = sqrt(linalg.det(covariance))
                bandwidth=KPDF.MPDFOptimumBandwidth(rr)

                if (self.kernel_estimator=="epanechnikov"):
                    lam = ones((shape(gdata)[0]),"float")
                    pdf2d = KPDF.MPDFEpanechnikov(rr,rr,bandwidth,lam,
                        invcovariance,renorm_constant)
                    lam = KPDF.MPDFAdapt(pdf2d)
                    pdf2d = KPDF.MPDFEpanechnikov(rr,gdata,bandwidth,lam,
                        invcovariance,renorm_constant)
                    pdf2d = reshape(pdf2d,(shape(kx)[0],shape(ky)[0]))
                else:
                    lam = ones((shape(gdata)[0]),"float")
                    pdf2d = KPDF.MPDFGaussian(rr,rr,bandwidth,lam,
                        invcovariance,renorm_constant)
                    lam = KPDF.MPDFAdapt(pdf2d)
                    pdf2d = KPDF.MPDFGaussian(rr,gdata,bandwidth,lam,
                        invcovariance,renorm_constant)
                    pdf2d = reshape(pdf2d,(shape(kx)[0],shape(ky)[0]))

            if self.verbose > 1:
                print "pdf2d=",pdf2d

            if sum(pdf2d.ravel()) > 0.:
                print "sum(pdf2d.ravel())",sum(pdf2d.ravel())
                karray = normalize(pdf2d,1.0)
                print "sum(karray.ravel())", sum(karray.ravel())
            else:
                karray = zeros(shape(pdf2d))

            # Normalize the kernel (accounting for incompleteness)
            if nmodel_bin <= 4: # if number of points used in kernel is less than 4
                completeness = -1.
            karray = karray * completeness
            if self.verbose:
                for i in range(len(bx0_bin)):
                    print "%7.3g %7.3g |" % (bx0_bin[i],bx1_bin[i]),
                print "%7d %7.3f " % (nmodel_bin,completeness)
                print '\n'

            kclass = kernel(bx0_bin,bx1_bin,self.pixdx,completeness,karray)
            kclass.completeness = completeness
            kclass.ninput_bin = bin_npoints
            kclass.nmodel_bin = nmodel_bin # number of points recovered by SE in the bin
            kclass.ngood_bin = n_goodfit_gf_bin
            kclass.covariance = covariance
            kclass.lam = lam
            kgrid.kernels[tuple(bxi)] = kclass        # Add the kernel to the grid

        print "len(kgrid.kernels)", len(kgrid.kernels) 
        self.kgrid = kgrid
        #return kgrid

    def write(self,filename):
        #self2 = copy.deepcopy(self)
        f = open(filename,'wb')
        cPickle.dump(self,f,2)
        f.close()


        
