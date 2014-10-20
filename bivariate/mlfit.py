#!/usr/bin/env python
# 
# Version 1.0.0 H. Ferguson 9/15/03
# Version 1.1.0 H. Ferguson 9/22/03 -- Added draw_from_pdf()
# Version 1.2.0 H. Ferguson 9/24/03 -- Revised to add verbose flags and 
#                                      to pass simplex fitting arguments
# Version 1.3.0 H. Ferguson 9/24/03 -- Added Poisson likelihood options, floor
# Version 2.0.0 H. Ferguson 1/6/05  -- Added routines to fit multiple data
#                                      sets
# Version 2.1.0 K-H Huang   5/21/09 -- Modify to use numpy instead of numarray

# Last Updated: 07/06/09  K-H Huang
"""Maximum-Likelihood Fit Routines for fitting a finely gridded 
multidimensional model to unbinned data. 
"""

__author__ =  'Henry C. Ferguson'
__version__=  '2.0.0'


from numpy import *
from scipy import *
from mlutil import *
from pygoods import *
import sys
from scipy.optimize import anneal
from uniparam import datalimits,pixdx
import bivariate_lf as bl


def loglikelihood(data,inputmodel,datalimits,floor=1.e-50):
    """Return the -sum(log(likelhood))

       Arguments:
       data - data array of shape (ndata,ndim)
       inputmodel - 2-D model array (pdf normalized to same total number) 
    """
    if shape(data)[0]!=2:
        raise ValueError, "input data array has incorrect shape"
    ndata = shape(data)[1] # number of data points
    model = maximum(inputmodel.model,floor)
    datagrid = grid_data(data,datalimits,inputmodel.pixdx)
    logl = sum((datagrid*log(model)).ravel())
    return -logl

#def mlfunc(parameters,data,func,func_otherarg,func_kwarg,
#            datalimits,floor,verbose):
#def mlfunc(parameters,alpha,data,func,func_otherarg,func_kwarg,
#           floor=1.e-50,verbose=2):
def mlfunc(parameters,data,restlimits,restpixdx,modellimits,obspixdx,datalimits,
           kgrid,norm,dropdic,kcorrfile,fft,meankcorr,verbose=1,steps="C"):
    """ Compute the loglikelihood of a data given a model. 

        Arguments:
        parameters -- parameters of the model
        #data -- masked data array (i.e. already gridded)
        ndata -- number of data points
        func_arg -- list argument of func
        func_kwarg -- keyword argument of func
        limits -- list of upper & lower limits
        nbins -- list of number of bins for each dimension
    """
    #if verbose:
    #    print parameters

    #pp = concatenate(([alpha],parameters))
    # If func == bl.bivariate_lf:
    #   - pp: parameters (alpha,Mstar,logr0,sigma,beta)
    #   - func_otherarg[0]=restlimits: restframe model limits
    #   - func_otherarg[1]=restpixdx: restframe model pixel size
    #   - func_otherarg[2]=modellimits: observed model limits (only cover the data points)
    #   - func_otherarg[3]=obspixdx: observed model pixel size
    #   - func_kwarg[0]=kgrid: transfer function kernel grid 
    #   - func_kwarg[1]=dropdic: dropout distribution 
    
    model = bl.bivariate_lf(parameters,restlimits,restpixdx,modellimits,obspixdx,datalimits,
            kgrid=kgrid,fft=fft,dropdic=dropdic,kcorrfile=kcorrfile,norm=norm,
            meankcorr=meankcorr,M0=parameters[1],steps=steps)
    #model = func(pp,*func_otherarg,**func_kwarg)
    model_array = maximum(model.model,1.e-50)     # get rid of negative values
    #pdf = model_array.copy()
    normfactor = 1.0/sum(model_array.ravel())

    # normalize the model within datalimits to a total of 1.0
    model.model = model_array * normfactor    # Normalize

    # Calculates log(likelihood)
    logl =  loglikelihood(data,model,datalimits) 

    if verbose:
        print "pars, log(l): ", parameters, logl
    return logl

def printlogl(parameters,data,restlimits,obslimits,datalimits,pixdx,kgrid=None,
              dropdic='bdrops.p',kcorrfile='kcorr1700_i.p'):
    #print shape(data)
    model = bl.bivariate_lf(parameters,restlimits,pixdx,obslimits,pixdx,kgrid=kgrid,
            dropdic=dropdic,kcorrfile=kcorrfile)
    idmin = findindex(model.limits,model.pixdx,datalimits[:,0])[0]
    idmax = findindex(model.limits,model.pixdx,datalimits[:,1])[0]
    #print idmin,idmax
    pdf = model.model[idmin[0]:idmax[0],idmin[1]:idmax[1]]
    normfactor = 1.0/sum(pdf.ravel())
    print normfactor
    model.model = maximum(model.model,1.e-50) * normfactor
    logl = zeros(shape(data)[1],'float')
    for i in range(shape(data)[1]):
        pts = data[:,i]
        #id = floor((pts-obslimits[:,0])/0.02).astype('int')
        id = findindex(obslimits,pixdx,pts)[0]
        logl[i] = -1.*log(model.model[id[0],id[1]])
    logl_tot = sum(logl)
    print "total logl", logl_tot
    for i in range(shape(data)[1]):
        print "%10s  %.2f  %.2f" % (data[:,i],logl[i],(logl[i]/logl_tot*100.))

### ===== Functions below are not used anymore ================= ###

#----------------------------------------------------------------------------
# Fitting Routines
#----------------------------------------------------------------------------

# Simplex 
# This routine probably will not fit all applications, but is 
# general enough that it makes sense to include it here. 
# 
simplex_defaults = {
    "xtol":1.e-6,
    "ftol":1.e-6,
    "maxiter":1000,
    "full_output":1
}    
def mlfitsimplex(datapoints,func,limits,nbins,guess,
    simplexargs=simplex_defaults,verbose=0,floor=1.e-50):
    """ Example maximum likelihood fit. Specific options for optimize.fmin
        should be adjusted for a specific case.

        Arguments:
        datapoints  -- an array of data of shape d,n where d is the number 
                of dimensions and n is the number of data points.
        func -- function which returns an array of model densities
                gridded as described by limits,nbins
        limits -- the limits of the coordinates
	nbins -- the number of bins in each axis
        guess -- the initial values of the free parameters.
        simplexargs -- options to pass to scipy simplex routine
        verbose -- 1 == print some diagnostics, 0 = be silent
    """

    # Bin the data
    ndata_in = shape(datapoints)[-1:][0]
    (data,ndata) = grid_data(datapoints,limits,nbins,verbose=verbose)     
    if verbose:
        print "Fitting %d data points (%d masked)" % (ndata,ndata_in-ndata)
    args = (data,ndata,func,limits,nbins,floor)
    simplexargs["args"] = args	  # Add to dictionary for passing keyword arguments
    simplexargs["disp"] = verbose
    
    # Do the fit
    (r,fval,iter,funcalls,warnflag) = apply(
            optimize.fmin,(mlfunc,(guess)),simplexargs)

    # Print the results
    if verbose:
        print "(r,fval,iter,funcalls,warnflag)"
        print r,fval,iter,funcalls,warnflag
    return (r,fval,iter,funcalls,warnflag)

# Annealing 
# This routine probably will not fit all applications, but is 
# general enough that it makes sense to include it here. 
# 
def mlfitanneal(datapoints,func,limits,nbins,guess,verbose=0,annealargs={}):
    """ Example maximum likelihood fit. Specific options for optimize.anneal
        should be adjusted for a specific case.

        Arguments:
        datapoints  -- an array of data of shape d,n where d is the number 
                of dimensions and n is the number of data points.
        func -- function which returns an array of model densities
                gridded as described by limits,nbins
        limits -- the limits of the coordinates
	nbins -- the number of bins in each axis
        guess -- the initial values of the free parameters.
    """

    # Bin the data
    ndata_in = shape(datapoints)[-1:][0]
    (data,ndata) = grid_data(datapoints,limits,nbins,verbose=verbose)     
    if verbose:
        print "Fitting %d data points (%d masked)" % (ndata,ndata_in-ndata)
    args = (data,ndata,func,limits,nbins)
    annealargs["args"] = args

    # Do the fit
    (r,jmin,T,feval,iter,accept,retval) = apply(
          anneal.anneal,(mlfunc,(guess)),annealargs)

    # Print the results
    if verbose:
        print "(r,jmin,T,feval,iter,accept,retval)"
        print r,jmin,T,feval,iter,accept,retval

    return (r,jmin,T,feval,iter,accept,retval)

#------------------------------------------------------------------------
# Markov chains
#------------------------------------------------------------------------
def mcmc_gauss(datapoints,func,limits,nbins,bestfit,sigma,
         nsamp,nburn=1000,verbose=0,floor=1.e-50):
    """ Metropolis Monte-Carlo Markov chain on a Gaussian proposal distribution 
        with given mean & sigma in each axis
        
        Arguments:
        datapoints  -- an array of data of shape d,n where d is the number
                of dimensions and n is the number of data points.
        func -- function which returns an array of model densities
                gridded as described by limits,nbins
        limits -- the limits of the coordinates
        nbins -- the number of bins in each axis
        bestfit -- the estimated best-fit values of the free parameters.
        sigma -- the proposal distribution sigma values for the free parameters.
        nsamp -- how many samples to take after burn in
        nburn -- how many samples to use for the burn in
    """           
    
    # Initialize
    (data,ndata) = grid_data(datapoints,limits,nbins,verbose=verbose)
    parameters = bestfit.copy()
    npar = shape(bestfit)[0]
    for i in range(len(bestfit)):
        parameters[i] = random.gauss(bestfit[i],sigma[i])
    mlogl = mlfunc(parameters,data,ndata,func,limits,nbins,floor=floor)
    logl = -mlogl
    burnresults = zeros((nburn,npar+1))*0.
    sampleresults = zeros((nsamp,npar+1))*0.
    currentpar = parameters.copy()
    currentresults = zeros(npar+1)*0.
  
    niter = 0
    while niter < nburn:
        for i in range(len(bestfit)):
            parameters[i] = random.gauss(currentpar[i],sigma[i])
        mlogl = mlfunc(parameters,data,ndata,func,limits,nbins,floor=floor)
        newlogl = -mlogl
        currentresults[:npar] = currentpar
        currentresults[npar] = newlogl
        accept = 0
        if newlogl >= logl:
            accept=1  
        if newlogl < logl:
            x = random.random()
            if newlogl-logl < -100:
                p = 0.
            else:
                p = exp(newlogl-logl)
            if x < p:
                accept = 1
        if accept:
            print "burn ", parameters, accept, niter, logl 
            currentpar = parameters.copy()
            burnresults[niter] = currentresults
            logl = newlogl
            niter += 1   

    niter = 0
    while niter < nsamp:
        for i in range(len(bestfit)):
            parameters[i] = random.gauss(currentpar[i],sigma[i])
        mlogl = mlfunc(parameters,data,ndata,func,limits,nbins,floor=floor)
        newlogl = -mlogl
        currentresults[:npar] = currentpar
        currentresults[npar] = newlogl
        accept = 0
        if newlogl >= logl:
            accept=1  
        if newlogl < logl:
            x = random.random()
            if newlogl-logl < -100:
                p = 0.
            else:
                p = exp(newlogl-logl)
            if x < p:
                accept = 1
        if accept:
            print "sample ",parameters, accept, niter, logl 
            currentpar = parameters.copy()
            sampleresults[niter] = currentresults
            logl = newlogl
            niter += 1   

    return (burnresults,sampleresults)

def mcmc_gauss_lim(datapoints,func,limits,nbins,bestfit,sigma,
         nsamp,nburn=1000,lower=[],upper=[], 
         verbose=0,floor=1.e-50,burnfile="",sampfile="",
         output_format=""):
    """ Metropolis Monte-Carlo Markov chain on a Gaussian proposal distribution 
        with given mean & sigma in each axis
        
        Arguments:
        datapoints  -- an array of data of shape d,n where d is the number
                of dimensions and n is the number of data points.
        func -- function which returns an array of model densities
                gridded as described by limits,nbins
        limits -- the limits of the coordinates
        nbins -- the number of bins in each axis
        bestfit -- the estimated best-fit values of the free parameters.
        sigma -- the proposal distribution sigma values for the free parameters.
        nsamp -- how many samples to take after burn in
        nburn -- how many samples to use for the burn in
        burnfile -- optional outputfile for burn in 
        sampfile -- optional outputfile for samples
        ouptut_format -- format of parameters if output is desired
    """           
    
   
    # Initialize
    if len(output_format)>0:
       print "lower:   "+output_format % tuple(lower)
       print "initial: "+output_format % tuple(bestfit)
       print "upper:   "+output_format % tuple(upper)
       print "sigma:   "+output_format % tuple(sigma)
    else:
       print "lower:   ",lower
       print "initial: ",bestfit
       print "upper:   ",upper
       print "sigma:   ",sigma


    if len(burnfile) > 0:
       print "burnfile = ", burnfile
       fburn = open(burnfile,'w')
    if len(sampfile) > 0:
       print "sampfile =", sampfile
       fsamp = open(sampfile,'w')
    if len(lower) > 0:
        apply_limits = 1
    else:
        apply_limits = 0
    (data,ndata) = grid_data(datapoints,limits,nbins,verbose=verbose)
    parameters = bestfit.copy()
    npar = shape(bestfit)[0]
    print "apply_limits = ",apply_limits
    if apply_limits:
        for i in range(len(bestfit)):
            parameters[i] = lower[i]-1.
            while (parameters[i] < lower[i] or parameters[i] > upper[i]):
                 parameters[i] = random.gauss(bestfit[i],sigma[i])
    else:
        for i in range(len(bestfit)):
            parameters[i] = random.gauss(bestfit[i],sigma[i])
    mlogl = mlfunc(parameters,data,ndata,func,limits,nbins,floor=floor)
    logl = -mlogl
    burnresults = zeros((nburn,npar+1))*0.
    sampleresults = zeros((nsamp,npar+1))*0.
    currentpar = parameters.copy()
    currentresults = zeros(npar+1)*0.

    niter = 0
    while niter < nburn:
        if apply_limits:
            for i in range(len(bestfit)):
                parameters[i] = lower[i]-1.
                while (parameters[i] < lower[i] or parameters[i] > upper[i]):
                    parameters[i] = random.gauss(currentpar[i],sigma[i])
        else:
            for i in range(len(bestfit)):
                parameters[i] = random.gauss(currentpar[i],sigma[i])

        mlogl = mlfunc(parameters,data,ndata,func,limits,nbins,floor=floor)
        newlogl = -mlogl
        currentresults[:npar] = currentpar
        currentresults[npar] = newlogl
        accept = 0
        if newlogl >= logl:
            accept=1  
        if newlogl < logl:
            x = random.random()
            if newlogl-logl < -100:
                p = 0.
            else:
                p = exp(newlogl-logl)
            if x < p:
                accept = 1
        if accept:
            if len(output_format) > 0:
                line = "%5d %1d %10.2f " % (niter,accept,logl)
                line = line + output_format % tuple(array(parameters).tolist())
                print "burn %s" % line
            else:
                print "burn ", parameters, accept, niter, logl 
            if len(burnfile) > 0:
                line = line + "\n"
                fburn.write(line)
                fburn.flush()
            currentpar = parameters.copy()
            burnresults[niter] = currentresults
            logl = newlogl
            niter += 1   
    fburn.close()

    niter = 0
    while niter < nsamp:
        if apply_limits:
            for i in range(len(bestfit)):
                parameters[i] = lower[i]-1.
                while (parameters[i] < lower[i] or parameters[i] > upper[i]):
                    parameters[i] = random.gauss(currentpar[i],sigma[i])
        else:
            for i in range(len(bestfit)):
                parameters[i] = random.gauss(currentpar[i],sigma[i])
        mlogl = mlfunc(parameters,data,ndata,func,limits,nbins,floor=floor)
        newlogl = -mlogl
        currentresults[:npar] = currentpar
        currentresults[npar] = newlogl
        accept = 0
        if newlogl >= logl:
            accept=1  
        if newlogl < logl:
            x = random.random()
            if newlogl-logl < -100:
                p = 0.
            else:
                p = exp(newlogl-logl)
            if x < p:
                accept = 1
        if accept:
            if len(output_format) > 0:
                line = "%5d %1d %10.2f " % (niter,accept,logl)
                line = line + output_format % tuple(array(parameters).tolist())
                print "samp %s" % line
            else:
                print "sample ", parameters, accept, niter, logl 
            if len(sampfile) > 0:
                line = line + "\n"
                fsamp.write(line)
                fsamp.flush()
            currentpar = parameters.copy()
            sampleresults[niter] = currentresults
            logl = newlogl
            niter += 1   
    fsamp.close()

    return (burnresults,sampleresults)

#------------------------------------------------------------------------
# Simultaneous fits of multiple distributions
#------------------------------------------------------------------------

# The only thing tricky here is the normalization. 
def mlfunc_multi(parameters,data_dict,ndata_dict,func,
                 limits_dict,nbins_dict,relnorm_dict,floor=1.e-50):
    """ Compute the loglikelihood of a data given a model. 
        The log likelihood is the sum of the log likelihoods of
        the individual data sets compared to the model.

        Arguments:
        The _dict arguments are dictionaries of items, one 
        for each data set
        parameters -- parameters of the model
        data -- masked data array (i.e. already gridded)
        ndata -- number of data points
        func -- routine to call to compute model from parameters
                This function returns the model and *two* normalization factors
        limits -- list of upper & lower limits
        nbins -- list of number of bins for each dimension
        relnorm -- normalization factor to put all models to same density
        floor -- minimum density allowed for model 
    """
    model_dict = {}
    completeness = {}
    eta = {}
    final_norm = {}
    zero = data_dict.keys()[0]
    
    # Get the models & normalization factors
    print "parameters: ",parameters 
    for n in data_dict.keys():
       ndata = ndata_dict[n]
       limits = limits_dict[n]
       nbins = nbins_dict[n]
       # Compute the model and normalization factors without & with 
       # the transfer function
       model_dict[n],norm,completeness[n] = apply(func,
                                                 (parameters,limits,nbins,n))
       integral = 1./norm
       eta[n] = integral
       print "%d: sum(model), relnorm[n], eta[n] " % n,sum(model_dict[n].ravel()),relnorm_dict[n], eta[n] 
    print "Completeness: ",completeness
    
    # First, renormalize based on ratio of predicted counts to observed
    # counts in the *first sample*. 
    ndata0 = ndata_dict[zero]
    b = ndata0/(completeness[zero]*eta[zero])
    
    total_pred = 0.
    total_obs = 0.
    pcounts = {}
    for n in data_dict.keys():
       pcounts[n] = b*relnorm_dict[n]*eta[n]*completeness[n]
       total_pred += pcounts[n]
       total_obs += ndata_dict[n]

    for n in data_dict.keys():
        final_norm[n] = (total_obs/total_pred)*b*relnorm_dict[n]

    # Check the normalization
    print "Normalization:"  
    print "Dataset    Observed    Predicted1   Predicted2 final_norm"
    for n in data_dict.keys():
       predicted_counts = final_norm[n]*eta[n]*completeness[n]
       print "%5d " % (n),
       print "  %10.2f    %10.2f  %10.2f  %10.4f" % \
           (ndata_dict[n],pcounts[n],predicted_counts,final_norm[n])

    # Compute the log likelihoods and sum them
    logl = 0.
    total_pred = 0.
    for n in data_dict.keys():
        data = data_dict[n]
        model = final_norm[n]*model_dict[n]
        s = sum(model.ravel())
        total_pred += s
        print"n, sum(model), ndata[n] ", n,s,ndata_dict[n]
        logl += loglikelihood(data,model,floor=floor) 
    print "tot_obs, sum(models) ",total_obs, total_pred
    print "pars, log(l): ", parameters, logl
    return logl

def simplex_multi(datapoints_dict,func,limits_dict,nbins_dict,guess,
    relnorm_dict,
    simplexargs=simplex_defaults,verbose=0,floor=1.e-50):
    """ Example maximum likelihood fit. Specific options for optimize.fmin
        should be adjusted for a specific case.

        Arguments:
        datapoints_dict  -- Dictionary of arrays of data of shape d,n 
                where d is the number of dimensions and n is the 
                number of data points.
        func -- function which returns an array of model densities
                gridded as described by limits,nbins
        limits_dict -- the limits of the coordinates
	nbins_dict -- the number of bins in each axis
        guess -- the initial values of the free parameters.
        simplexargs -- options to pass to scipy simplex routine
        verbose -- 1 == print some diagnostics, 0 = be silent
    """

    # Bin the data
    #data_dict = {}
    #ndata_dict = {}
    #for n in datapoints_dict.keys():
    #   limits = limits_dict[n]
    #   nbins = nbins_dict[n]
    #   datapoints = datapoints_dict[n]
    #   ndata_in = shape(datapoints)[-1:][0]
    #   (data,ndata) = grid_data(datapoints,limits,nbins,verbose=verbose) 
    #   if verbose:
    #       print "Dataset %d Fitting %d data points (%d masked)" % (
    #              n,ndata,ndata_in-ndata)
    #   data_dict[n] = data
    #   ndata_dict[n] = ndata
    #args = (data_dict,ndata_dict,func,limits_dict,nbins_dict,relnorm_dict,floor)
    #simplexargs["args"] = args  # Add to dictionary for passing keyword arguments
    #simplexargs["disp"] = verbose

    print "simplexargs: ",simplexargs
    
    # Do the fit
    #(r,fval,iter,funcalls,warnflag) = apply(
    #        optimize.fmin,(mlfunc_multi,(guess)),simplexargs)
    (r,fval,iter,funcalls,warnflag) = optimize.fmin(func,guess,**simplexargs)
    
    # Print the results
    if verbose:
        print "(r,fval,iter,funcalls,warnflag)"
        print r,fval,iter,funcalls,warnflag
    return (r,fval,iter,funcalls,warnflag)

