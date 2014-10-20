#!/usr/bin/env python

import numpy as np

def draw_from_pdf(ntot, pdf_array, x0_array):
    """
    Draws ntot number between x0_array[0], x0_array[-1]+dx, with the specified
    probability distribution pdf_array.
    Assumes that x0_array is a regular grid.
    """
    dx = x0_array[1] - x0_array[0]
    # Build a cumulative PDF
    pdf_cum = np.cumsum(pdf_array)
    # Make the last element = 1.0
    pdf_cum = pdf_cum / pdf_cum[-1] 
    # Generate an array of random numbers between 0 and 1
    randx = np.random.uniform(0., 1., size=ntot)
    # Search the indices of these random numbers in the cumulative PDF
    irandx = np.searchsorted(pdf_cum, randx)
    # force minimum index to be 0
    irandx = np.maximum(0, irandx)
    irandx = np.minimum(len(pdf_cum)-1, irandx)
    # Get the values of x at these indices
    xrand = x0_array.take(irandx)
    # Perturb these values a little bit
    dx_array = np.random.uniform(0., dx, size=ntot)
    xrand = xrand + dx_array
    return xrand

def draw_from_pdf_func_old(ntot,pdf,limits,ndim=1,verbose=0):
    """ Return ntot data points drawn from a model distribution function.

        Arguments:
        ntot -- number of data values to create
        pdf -- a FUNCTION that returns the pdf P(x)
               pdf(x) returns the value of PDF at x
        limits -- the limits of drawn values
    """
    if ndim == 1:
        datapoints = np.zeros(ntot, 'float')
    else:
        datapoints = np.zeros((ndim,ntot),"float")
    index = np.zeros(ndim,"int")
    # normalize the pdf
    n = 0
    for n in range(ntot):                # Populate with random draws
        drawn = 0
        while drawn == 0:
            #ilow = np.ones(ndimensions) * -0.5
            #ihigh = np.array(np.shape(pdf)) - 0.51
            #index = np.around(np.random.uniform(ilow, ihigh)).astype('int')
            #index = np.around(np.random.uniform([-0.5,-0.5],np.array(np.shape(pdf))-0.51)).astype('int')
            x = np.random.uniform(limits[0], limits[1])
            y = np.random.random()  
            ydraw = pdf(x)
            if ydraw > y:  # then draw a point from this bin
                drawn = 1
                if ndim == 1:
                    datapoints[n] = x
                else:
                    datapoints[:,n] = x
                
    return datapoints