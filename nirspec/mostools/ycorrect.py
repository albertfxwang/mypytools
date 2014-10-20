"""
Module to determine the y-distortion of multi-slit mask images. This code
  hasn't been looked at in a *very* long time....
"""

import scipy,special_functions
from scipy import ndimage,signal

def ycorrect(data,dostar=True,order=4):
    """
    ycorrect(data)

    Inputs:
      data - a flatfield image of the mask

    Outputs:
      true_coeffs - A polynomial describing the transformation:
                    y_straight = f(x_ccd,y_ccd)
      map_coeffs  - A polynomial describing the transformation:
                    y_ccd = f(x_cdd,y_straight)
    """

    # Parameters
    SUMWIDTH = 41    # Width of summing over columns

    y_axis = data.shape[0]
    x_axis = data.shape[1]


    central = x_axis/2

    x_min_orig = central - SUMWIDTH/2
    x_max_orig = central + SUMWIDTH/2

    # Find the 'holes' in the center of the mask to use as the reference
    #   position.
    midcol = data[:,x_min_orig:x_max_orig].mean(axis=1)
    central_edges,threshold,star_cutoff = find_holes(midcol,dostar)

    # transform_table would be easier to use as a list....
    transform_table = scipy.zeros((1,3),'f4')
    index = 0
    for peak in central_edges:
        if index:
            transform_table.resize((index+1,3))
        transform_table[index,0] = central
        transform_table[index,1] = peak
        transform_table[index,2] = peak

        index += 1

    offset = scipy.zeros(len(central_edges))

    x_min = x_min_orig
    x_max = x_max_orig
    current_column = central

    while current_column>SUMWIDTH + 20:
        current_column = current_column - SUMWIDTH - 10
        x_min = x_min - SUMWIDTH - 10
        x_max = x_max - SUMWIDTH - 10

        comp_array = data[:,x_min:x_max].mean(axis=1)
        comp_array.clip(min=-1000.,max=star_cutoff)
        derivative = deriv_1d(comp_array)
        derivative = ndimage.gaussian_filter1d(derivative,3)
        derivative = abs(derivative)

        for i in range(offset.size):
            if scipy.isnan(offset[i]):
                continue
            ref = central_edges[i] + offset[i]

            start = int(ref) - 6
            end = start + 13

            if derivative[start:end].max()<threshold:
                offset[i] = scipy.nan
                continue

            fit = find_peak(derivative[start:end])

            # If the fit has crazy parameters, skip it
            if(fit[2]<0 or fit[2]>13 or fit[3]<1 or fit[3]>6):
                offset[i] = scipy.nan
                continue

            peak = fit[2]+float(start)            
            offset[i] = peak - central_edges[i]

            transform_table.resize((index+1,3))

            transform_table[index,0] = current_column
            transform_table[index,1] = central_edges[i]
            transform_table[index,2] = peak

            index += 1

    offset = scipy.zeros(offset.size)

    x_min = x_min_orig
    x_max = x_max_orig
    current_column = central
    while current_column<x_axis - SUMWIDTH - 19:
        current_column = current_column + SUMWIDTH + 10
        x_min = x_min + SUMWIDTH + 10
        x_max = x_max + SUMWIDTH + 10

        comp_array = data[:,x_min:x_max].mean(axis=1)
        comp_array.clip(min=-1000.,max=star_cutoff)
        derivative = deriv_1d(comp_array)
        derivative = ndimage.gaussian_filter1d(derivative,3)
        derivative = abs(derivative)

        for i in range(offset.size):
            if scipy.isnan(offset[i]):
                continue
            ref = central_edges[i] + offset[i]
 
            start = int(round(ref)) - 6
            end = start + 13
 
            if derivative[start:end].max()<threshold:
                offset[i] = scipy.nan
                continue
 
            fit = find_peak(derivative[start:end])

            if(fit[2]<0 or fit[2]>13 or fit[3]<1 or fit[3]>6):
                offset[i] = scipy.nan
                continue
 
            peak = fit[2]+float(start)
            offset[i] = peak - central_edges[i]

            transform_table.resize((index+1,3))

            transform_table[index,0] = current_column
            transform_table[index,1] = central_edges[i]
            transform_table[index,2] = peak

            index += 1

    true_coeffs = special_functions.lsqfit(transform_table,"chebyshev",order,order)

    temp = transform_table[:,1].copy()
    transform_table[:,1] = transform_table[:,2].copy()
    transform_table[:,2] = temp.copy()

    map_coeffs = special_functions.lsqfit(transform_table,"chebyshev",order,order)

    return true_coeffs,map_coeffs


"""
Supporting Functions
"""
def find_holes(data,dostar=True):
    """
    find_holes(data)

    Finds the edges of slits to use as references for determining the
      distortion solution. I can't remember why it's named this, and I'm
      sure I can come up with a more elegant solution than what's here, but
      it seems to work, so I'll leave it.

    Input:
      data - central column of flatfield

    Output:
      edges       - list containing the locations of the slit edges
      threshold   - the threshold used to select valid edges
      star_cutoff - the plateau of valid data
    """
    sample = data.copy()
    size = sample.size

    # Here's a little hack to "flatten" star boxes
    tmp = scipy.sort(sample)
    star_cutoff = scipy.median(tmp[-30:-10])*0.6
    if dostar is True:
        sample[sample>star_cutoff] = star_cutoff

    derivative = deriv_1d(sample)
    derivative = ndimage.gaussian_filter1d(derivative,3)
    derivative = abs(derivative)

    def clip(arr,nsig=3.):
        a = arr.flatten()
        a.sort()
        a = a[a.size*0.05:a.size*0.75]
        m,s,l = a.mean(),a.std(),a.size
        while 1:
            a = a[abs(a-m)<s*nsig]
            if a.size==l:
                return m,s
            m,s,l = a.mean(),a.std(),a.size

    tmp = scipy.sort(derivative)
    avg = scipy.median(tmp[size/8:size*3/8])
    sigma = tmp[size/8:size*3/8].std()
#    avg,sigma = clip(derivative)

    threshold = avg + sigma*100.

    edge = []

    count = 0
    while derivative.max()>threshold:
        start = derivative.argmax()-7
        end = derivative.argmax()+8

        if start<0:
            start = 0
        if end>derivative.size:
            end = derivative.size

        fit = find_peak(derivative[start:end])

        if start>7 and end<derivative.size-7:
            edge.append(float(start) + fit[2])

        start -= 3
        end += 3

        if start<0:
            start = 0
        if end>derivative.size:
            end = derivative.size

        derivative[start:end] = 0.

    edge.sort()
    return edge,threshold,star_cutoff


"""
A simple weighted average would probably work too....
"""
def find_peak(data):
    size = data.shape[0]
    fit_data = scipy.reshape(scipy.arange(0,size,0.5),(size,2))
    fit_data[:,1] = data.copy()

    fit = scipy.zeros(4)
    fit[0] = 0.
    fit[1] = data.max()
    fit[2] = data.argmax()
    fit[3] = 2.
    fit,val = special_functions.ngaussfit(fit_data,fit)

    return fit

def deriv_1d(data):
    out = signal.convolve(data,[-0.5,0.,0.5],mode='same')
    size = data.size

    # Deal with the ends properly
    out[0] = data[1] - data[0]
    out[size-1] = data[size-1]-data[size-2]

    return out
