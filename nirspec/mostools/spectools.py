"""
A suite of functions for manipulating 2d and 1d spectra.
"""

import scipy,pyfits


def resample1d(data,coeffs,axis,const=0.0,offset=0.):
    """
    resample1d(data,coeffs,axis,const=0.0,offset=0.)

    Resamples 2d spectra in x or y direction.

    Inputs:
      data   - input data array
      coeffs - polynomial that describes output
      axis   - axis to resample ("x" or "y")
      const  - padding constant for the interpolation
      offset - amount by which to offset the coefficient solution

    Output:
      resampled data
    """
    from special_functions import genfunc
    coords = array_coords(data.shape)
    y = coords[0].flatten()
    x = coords[1].flatten()
    if axis=="Y" or axis=="y" or axis==1:
        y += offset
        coords[0] = genfunc(x,y,coeffs).reshape(data.shape)-offset
    else:
        x += offset
        coords[1] = genfunc(x,y,coeffs).reshape(data.shape)-offset
    del x,y
    return scipy.ndimage.map_coordinates(data,coords,cval=const,output=scipy.float64,order=5)

def array_coords(shape):
    """
    Faster version of scipy.indices()
    """
    return scipy.indices(shape).astype(scipy.float32)
    y = shape[0]
    x = shape[1]
    out = scipy.empty((2,y,x))
    t = scipy.arange(y,dtype='f8')
    out[0] = scipy.tile(t,(x,1)).T
    t = scipy.arange(x,dtype='f8')
    out[1] = scipy.tile(t,(y,1))
    return out

def resampley(data,ycoords,yoffset=0.,cval=0.,mode="constant",slice=None):
    """
    resampley(data,ycoords,yoffset=0.,cval=0.,mode="constant",slice=None)

    Resamples 2d spectra in the spatial direction to remove the grating
      smile. Allows a coordinate array or polynomial to be passed.

    Inputs:
      data    - data to be resampled
      ycoords - output y-coordinates or a polynomial to describe output
      yoffset - optional offset wrt the polynomial or coordinate array
      cval    - interpolation boundary constant
      mode    - method for interpolating boundaries
      slice   - sub-slice to resample to (ie if the output shape is not the
                   same as the input shape)
    """

    from scipy import ndimage

    try:
        coords = array_coords(ycoords.shape)
    except:
        from special_functions import genfunc
        coords = array_coords(data.shape)
        x = coords[1].flatten()
        y = coords[0].flatten()+yoffset
        ycoords = genfunc(x,y,ycoords).reshape(coords[0].shape)
        if slice is not None:
            ycoords = ycoords[slice].copy()
        coords = array_coords(ycoords.shape)
        del x,y
    coords[0] = ycoords.astype(scipy.float64)-yoffset
    return ndimage.map_coordinates(data,coords,output=scipy.float64,mode=mode,cval=cval,order=5)

def cutout(file,outname,slit_number,plane=0):
    """
    cutout(file,outname,slit_number,plane=0)

    Cut out a slit from a 2d mask image.

    Inputs:
      file        - input 2d mask image
      outname     - name of output FITS file containing slit
      slit_number - the number of the slit to cut out (with 1 as the bottom)
      plane       - the plane if file is a FITS cube
    """

    f = pyfits.open(file)
    data = f[0].data.copy()
    hdr = f[0].header.copy()
    data[scipy.isnan(data)] = 0.

    if len(data.shape)==3:
        data = data[plane]
    slit,start,bottom = cut_slit(data,slit_number)
    crval = hdr['CRPIX1']-start

    hdu = pyfits.PrimaryHDU(slit)
    hdu.header = hdr
    hdu.header.update('CRPIX1',crval)
    hdu.verify('fix')

    hdu.writeto(outname)
    

def get_slit(data,slit_number):
    """
    get_slit(data,slit_number)

    Return the slit as an array.

    Inputs:
      data        - 2d array of mask image
      slit_number - the number of the slit to cut out (with 1 as the bottom)

    Outputs:
       2d array of slit data
    """
    a,b,c = cut_slit(data,slit_number)
    return a

def cut_slit(data,slit_number):
    """
    cut_slit(data,slit_number)

    Find slit in 2d mask image.

    Inputs:
      data        - 2d array of mask image
          slit_number - the number of the slit to cut out (with 1 as the bottom)

    Outputs:
      2d array of slit data, left index of slit, bottom index of slit
    """
    data = data.copy()
    slit_number -= 1

    data[scipy.isnan(data)] = 0.

    slice = data.sum(axis=1)

    # Find bottom of good data
    bottom = 0
    while slice[bottom]==0:
        bottom += 1

    slice = scipy.trim_zeros(slice)

    zeros = scipy.where(slice==0)
    zeros = zeros[0]+bottom

    # Special case for images with only one slit
    if zeros.size==0:
        # Deal with upper zeros
        slice = data.sum(axis=1)
        top = bottom
        while top<slice.size and slice[top]!=0:
            top += 1

        slice = data.sum(axis=0)
        indx = scipy.where(slice!=0)
        start = indx[0][0]
        end = indx[0][-1]+1
        return data[bottom:top,start:end],start,bottom

        borders = []
        start = 0
        for i in range(zeros.size-1):
            if zeros[i]+1!=zeros[i+1]:
                borders.append(zeros[start:i+1].copy())
                start = i+1
        borders.append(zeros[start:].copy())
    if slit_number>len(borders):
        return scipy.asarray(0),0,0

    if slit_number==0:
        start = 0
        end = borders[0][0]
    elif slit_number==len(borders):
        start = borders[slit_number-1][4]+1
        end = slice.size
    else:
        start = borders[slit_number-1][4]+1
        end = borders[slit_number][0]

    data = data[start:end]
    bottom = start
    slice = data.sum(axis=0)
    indx = scipy.where(slice!=0)
    start = indx[0][0]
    end = indx[0][-1]+1
    return data[:,start:end],start,bottom


def logrebin(spectrum,crval,crpix,cd1):
    from scipy import interpolate
    # First determine the "best" pixel scale
    start = crval+cd1*(1-crpix)
    inwave = scipy.arange(start,start+spectrum.size*cd1,cd1)
    sampwave = scipy.log10(inwave)
    pixscale = (sampwave[-1]-sampwave[0])/sampwave.size
    outwave = scipy.arange(sampwave[0],sampwave[-1],pixscale)
    outwave = scipy.power(10,outwave)

    spline = interpolate.splrep(inwave,spectrum,s=0)
    newspec = interpolate.splev(outwave,spline)
    return outwave,newspec

def fits_logrebin(infile,outfile):
    import pyfits
    f = pyfits.open(infile)
    header = f[0].header
    data = f[0].data.copy()
    crval1 = header['crval1']
    crpix = header['crpix1']
    cd = header['cd1_1']
    outwave,outspec = rebin_log(data,crval1,crpix,cd)
    out = pyfits.PrimaryHDU(outspec)
    outwv = scipy.log10(outwave)
    start = outwv[0]
    delt = outwv[1]-outwv[0]
    out.header.update('CRVAL1',start)
    out.header.update('CRPIX1',1)
    out.header.update('CD1_1',delt)
    out.header.update('DC-FLAG',1)
    out.writeto(outfile)

def optimal_smooth(data,var,width=1.5):
    """
    optimal_smooth(data,var,width=1.5)

    Gaussian smoothing using variance information.

    Inputs:
      data  - 1d science spectrum
      var   - 1d variance spectrum
      width - width of Gaussian sigma in pixels

    Outputs:
      smoothed 1d science spectrum
    """
    from scipy import signal
    window = 8*width # Go out to 4sigma in the tails
    if window%2==0:
        window += 1
    w = signal.gaussian(window,width)
    tmp = data/var
    smooth = signal.convolve(tmp,w,'same')
    norm = signal.convolve(1./var,w,'same')
    return smooth/norm


def parse_hdr(header):
    """
    parse_hdr(header)

    Helper function that returns spectral WCS keywords from a FITS header.

    Inputs:
      header - FITS header

    Outputs:
      list containing crval, crpix, cd, and a flag indicating whether the
        wavelength is in lambda or log-lambda coordinates.
    """
    crval = header['crval1']
    try:
        crpix = header['crpix1']
    except:
        crpix = 1
    log = 0
    try:
        cd = header['cd1_1']
    except:
        cd = header['cdelt1']
    try:
        log = header['dc-flag']
    except:
        try:
            tmp = header['WFITTYPE']
            if tmp=='LOG_LINEAR':
                log = 1
        except:
            pass 
    return [crval,crpix,cd,log]

def wavelength(filename,ext=0):
    """
    wavelength(filename,ext=0)

    Creates an array representing the wavelength of a file based on the FITS
      WCS keywords.

    Inputs:
      filename - FITS file name
      ext      - extension number

    Outputs:
      1d array describing the wavelength of each pixel
    """
    f = pyfits.open(filename)
    hdr = f[ext].header

    hdr_info = parse_hdr(hdr)
    crval = hdr_info[0]
    crpix = hdr_info[1]
    cd = hdr_info[2]
    islog = hdr_info[3]
    npix = hdr['NAXIS1']

    start = crval+(1.-crpix)*cd
    wave = scipy.arange(start,start+npix*cd,cd)
    if islog:
        wave = scipy.power(10.,wave)
    return wave[0:npix]

def make_spec(spec,varspec,wave,outname,pos=-1,width=-1,clobber=False):
    """
    make_spec(spec,varspec,wave,outname)

    Create a FITS file describing an extracted spectrum.

    Inputs:
      spec    - science spectrum
      varspec - variance spectrum
      wave    - wavelength array describing wavelength of each science pixel
      outname - name of output FITS file

    Outputs:
      multi-extension FITS file with extension 1 containing the science
        spectrum, extension 2 has a smoothed version of the science
        spectrum, and extension 3 contains the variance spectrum.
    """

    smooth = optimal_smooth(spec,varspec,3)
    smooth[scipy.isnan(smooth)] = 0.

    hdu = pyfits.PrimaryHDU()
    hdu.header.update('CENTER',pos)
    hdu.header.update('WIDTH',width)
    hdulist = pyfits.HDUList([hdu])

    crval = wave[0]
    scale = wave[1]-wave[0]
    islog = False
    if scale!=wave[2]-wave[1]:
        test = scipy.log10(wave)
        if abs((test[1]-test[0])-(test[2]-test[1]))<1e-7:
            islog = True
            crval = test[0]
            scale = test[1]-test[0]
    for i in [spec,smooth,varspec]:
        thdu = pyfits.ImageHDU(i)
        thdu.header.update('CRVAL1',crval)
        thdu.header.update('CD1_1',scale)
        thdu.header.update('CRPIX1',1)
        thdu.header.update('CRVAL2',1)
        thdu.header.update('CD2_2',1)
        thdu.header.update('CRPIX2',1)
        if islog:
            thdu.header.update('CTYPE1','WAVE-LOG')
            thdu.header.update('DC-FLAG',1)
        else:
            thdu.header.update('CTYPE1','LINEAR')
        hdulist.append(thdu)

    hdulist.writeto(outname,clobber=clobber)

def coadd(infiles,oname,clobber=True):
    """
    coadd(infiles,oname)

    Combine spectra with inverse variance weighting
    """
    wave = wavelength(infiles[0],1)
    spec = wave*0.
    wht = wave*0.
    for file in infiles:
        t = pyfits.open(file)
        spec += t[1].data.copy()/t[3].data.copy()
        wht += 1./t[3].data.copy()
    spec = spec/wht
    wht = 1./wht

    make_spec(spec,wht,wave,oname,clobber=clobber)
