"""
Determine the spectral resolution.
"""

import scipy
from scipy import ndimage

import special_functions

def clipped(data,clip=3.5):
        d = data.copy()
        while 1:
                mean,std,size = d.mean(),d.std(),d.size
                d = d[abs(d-mean)<clip*std]
                if d.size==size or d.size==0:
                        return mean,std
        return mean,std


def get_resolution(file):
    import pyfits
    from mostools import spectools as st
    spec = pyfits.open(file)[3].data.copy()
    wave = st.wavelength(file,1)
    v = measure(spec,wave)
    if v==1:
        print 'No lines fit!'
        return 0
    return 299792./v[0]


def get_line_resolution(line,file):
    import pyfits
    from mostools import spectools as st
    spec = pyfits.open(file)[3].data.copy()
    wave = st.wavelength(file,1)
    data = spec.astype(scipy.float64)
    data[scipy.isnan(data)] = 0.
    size = data.size

    bg = ndimage.percentile_filter(data,20,55)
    bg = ndimage.gaussian_filter(bg,7)
    data -= bg

    pos = abs(wave-line).argmin()
    fitdata = scipy.zeros((15,2))
    fitdata[:,0] = wave[pos-7:pos+8]
    fitdata[:,1] = data[pos-7:pos+8]

    par = scipy.zeros(4)
    par[1] = data[pos]
    par[2] = line
    par[3] = 1.

    fit,chi2 = special_functions.ngaussfit(fitdata,par)

    pixsize = wave[pos]-wave[pos-1]
    if fit[3]>5.*pixsize or fit[3]<0.8*pixsize or abs(fit[2]-wave[pos])>1.5*fit[3]:
        print 'Could not fit for resolution of spectral line'
        print fit
        return 0

    v = fit[2]/fit[3]
    return 299792./v


def measure(spectrum,wave,nsig=30.,num=25,show=False):
    """
    measure(spectrum,num=25)

    Deterines the spectral resolution for a spectrum by measuring the width
      of arc or skylines.

    Inputs:
      spectrum - 1d spectrum
      wave     - array containing wavelengths of points in file
      num      - maximum number of lines to measure (not implemented)

    Outputs:
      median resolution of all lines extracted
    """

    data = spectrum.astype(scipy.float64)
    data[scipy.isnan(data)] = 0.
    size = data.size

    bg = ndimage.percentile_filter(data,20,35)
    bg = ndimage.gaussian_filter(bg,7)
    data -= bg

    avg,std = clipped(data,2.0)
    thresh = avg + nsig*std  # Only measure 30 sigma lines


    """ Identify isolated peaks. """
    mask = ndimage.maximum_filter(data,15)
    lines = scipy.where((mask==data)&(mask>thresh)&(mask<50000.))[0]
    count = 0

    vals = []
    locs = []
    for pos in lines:
        """ Don't use lines near the edges. """
        if pos<7 or pos+8>size:
            continue
        fitdata = scipy.zeros((15,2))
        fitdata[:,0] = wave[pos-7:pos+8]
        fitdata[:,1] = data[pos-7:pos+8]

        par = scipy.zeros(4)
        par[1] = data[pos]
        par[2] = wave[pos]
        par[3] = 1.

        fit,chi2 = special_functions.ngaussfit(fitdata,par)
        """
        Reject fits that were 'too' wide or narrow, or not near the
          expected center.
        """
        pixsize = wave[1]-wave[0]
        if fit[3]>5.*pixsize or fit[3]<0.8*pixsize or abs(fit[2]-wave[pos])>1.5*fit[3]:
            continue


        locs.append(fit[2])
        vals.append(fit[2]/fit[3])
        
        count += 1
    vals = scipy.asarray(vals)

    if vals.size==0:
        return 1.

    if show:
        import pylab
        pylab.plot(locs,vals)
        pylab.show()

    return scipy.median(vals),vals.size
