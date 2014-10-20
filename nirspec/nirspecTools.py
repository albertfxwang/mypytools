import pyfits,numpy
from scipy import ndimage,signal,interpolate,optimize
from mostools import spectools as st,specfit as spf
from mostools import ycorrect
import indexTricks as iT
import special_functions as sf
import nirspec,crsub
from math import cos,sin,pi


def nirspecOpen(filename, verbose=False):
    """
    python doesn't like the Keck headers....
    """
    if verbose:
        print "Reading file %s" % filename
    return pyfits.open(filename,ignore_missing_end=True)


def ycor(img, nstd=5.):
    """
    Follow star traces along the x direction to correct y-distortion
    """
    xvals = numpy.arange(img.shape[1])
    yvals = numpy.arange(img.shape[0])
    mid = spf.get_lines(numpy.arange(img.shape[0]),img[:,yvals.size/2],False,
                        nstd=nstd)
    ref = []
    for col in range(xvals.size):
        ref.append(spf.get_lines(yvals,img[:,col],False))
    soln = []
    for l in mid:
        off = 0.
        dy = []
        for col in xrange(xvals.size/2,xvals.max()):
            if len(ref[col])==0:
                continue
            p = ref[col] + off
            diff = abs(l-p)
            if diff.min()<1.5:
                pos = diff.argmin()
                dy.append([col,ref[col][pos]-l])
                off = l-ref[col][pos]
        off = 0.
        for col in xrange(xvals.size/2-1,-1,-1):
            if len(ref[col])==0:
                continue
            p = ref[col] + off
            diff = abs(l-p)
            if diff.min()<1.5:
                pos = diff.argmin()
                dy.append([col,ref[col][pos]-l])
                off = l-ref[col][pos]
        dy = numpy.asarray(dy)
        fit = sf.lsqfit(dy,'chebyshev',3)
        soln.append(sf.genfunc(xvals,0.,fit))
    soln = numpy.asarray(soln)
    soln = numpy.median(soln,0)
    y = numpy.indices(img.shape).astype(numpy.float32)[0]+soln
    return y,y-soln*2


def FFTFilter(img,width):
    """
    High pass filtering img
    """
    fft = numpy.fft.rfft(img)
    x = numpy.arange(fft.size).astype(numpy.float32)
    mod = numpy.fft.irfft(fft*x/(width+x))
    return mod


def clip(arr,nsig=3.5):
    a = arr.flatten()
    a.sort()
    a = a[a.size/20:a.size*0.95]
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<s*nsig]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


y,x = iT.coords((1200,1200))
c = cos(-84.9*pi/180.)
s = sin(-84.9*pi/180.)
X = (x-599.5)*c + (y-599.5)*s + 511.5
Y = (599.5-x)*s + (y-599.5)*c + 511.5
rotcoords = numpy.array([Y,X])
y,x = iT.coords((1024,1024))
c = cos(84.9*pi/180.)
s = sin(84.9*pi/180.)
X = (x-511.5)*c + (y-511.5)*s + 599.5
Y = (511.5-x)*s + (y-511.5)*c + 599.5
derotcoords = numpy.array([Y,X])
y,x = iT.coords((2400,2400))
c = cos(-84.9*pi/180.)
s = sin(-84.9*pi/180.)
X = (x-1199.5)*c + (y-1199.5)*s + 1023.5
Y = (1199.5-x)*s + (y-1199.5)*c + 1023.5
rotcoords2 = numpy.array([Y,X])
y,x = iT.coords((2048,2048))
c = cos(84.9*pi/180.)
s = sin(84.9*pi/180.)
X = (x-1023.5)*c + (y-1023.5)*s + 1199.5
Y = (1023.5-x)*s + (y-1023.5)*c + 1199.5
derotcoords2 = numpy.array([Y,X])
def rotate(img):
    return ndimage.map_coordinates(img,rotcoords)

def derotate(img):
    return ndimage.map_coordinates(img,derotcoords)

def rotate2(img):
    return ndimage.map_coordinates(img,rotcoords2)

def derotate2(img):
    return ndimage.map_coordinates(img,derotcoords2)


def illumination(flat,nsig=15,nhigh=200):
    s = clip(flat[500:600,300:400])[1]
    flatrot = rotate(flat)#ndimage.rotate(flat,84.9,order=5)
    tmp0 = numpy.where(flatrot>nsig*s,1,0)
    tmp1 = numpy.where(tmp0.sum(1)>nhigh)[0]
    low,high = tmp1.min(),tmp1.max()+1
    tmp2 = numpy.where(tmp0.sum(0)>nhigh/4)[0]
    left,right = tmp2.min(),tmp2.max()+1
    return low,high,left,right


def makeTraceIm(stars,nx,ny,save_output=False):
    #
    # Make trace image by stacking stars; if there are more than 2 then
    #   subtract the median
    #

    print ""
    print "Making trace image from star exposures"
    print "--------------------------------------"
    star = numpy.empty((len(stars),ny,nx))
    count = 0
    for img in stars:
        star[count] = nirspecOpen(img,verbose=True)[0].data.copy()
        count += 1
    del img
    bpm = numpy.median(star,0)

    star = star.sum(0)
    if count>2:
        star -= bpm*count

    if save_output:
        pyfits.PrimaryHDU(star).writeto('star.fits',clobber=True)

    return star

def getYsoln(star,nstars,low,high):
    # Find solution as function of y position if there are enough star traces
    print ""
    print "Finding y solution"
    star = star.repeat(2,0).repeat(2,1)
    star = rotate2(star)
    starTMP = star[low*2:high*2].copy()
    if nstars>3:
        ytrue2,ymap2 = ycorrect.ycorrect(starTMP,False,order=2)
        coords = numpy.indices(starTMP.shape).astype(numpy.float32)
        x = coords[1].flatten()
        y = coords[0].flatten()
        del coords
        yforw2 = sf.genfunc(x,y,ytrue2).reshape(starTMP.shape)
        yback2 = sf.genfunc(x,y,ymap2).reshape(starTMP.shape)
    # Not enough star traces; assume the distortion is constant over y
    else:
        yforw2,yback2 = ycor(starTMP,)
    # Repeat for non-oversampled case
    star = iT.resamp(star,2)[low:high]
    if nstars>3:
        ytrue,ymap = ycorrect.ycorrect(star,False,order=2)
        coords = numpy.indices(star.shape).astype(numpy.float32)
        x = coords[1].flatten()
        y = coords[0].flatten()
        del coords
        yforw = sf.genfunc(x,y,ytrue).reshape(star.shape)
        yback = sf.genfunc(x,y,ymap).reshape(star.shape)
    else:
        yforw,yback = ycor(star)
    return yforw,yback,yforw2,yback2#,ytrue2,ymap2


def getStraight(img,ysoln,low,resamp=True,mode='constant'):
    yforw = ysoln[2]
    tmp = img.repeat(2,0).repeat(2,1)
    x = iT.coords(yforw.shape)[1]
    c = numpy.array([yforw+low*2,x])
    X = ndimage.map_coordinates(rotcoords2[1],c)
    Y = ndimage.map_coordinates(rotcoords2[0],c)
    c = numpy.array([Y,X])
    tmp = ndimage.map_coordinates(tmp,c,mode=mode)
    if resamp==True:
         return iT.resamp(tmp,2)
    return tmp


def getFlatNorm(flat,ysoln,low,high):
    from scipy.stats import stats
    f = getStraight(flat,ysoln,low,False)
    ftmp = numpy.where(f==0,numpy.nan,f)
    norm = stats.nanmedian(ftmp[20:-20],0)
    norm = numpy.where((norm<=0)|numpy.isnan(norm),1.,norm)
    f /= norm
    img = numpy.ones((2400,2400))
    img[low*2:high*2] = st.resampley(f,ysoln[3],mode='nearest')
    return iT.resamp(derotate2(img),2)


def getWsoln(sky,x,wsolution,wmodel):
    def refine(p,x,d,sky,model):
        pars = {'coeff':numpy.atleast_2d(p).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,pars)
        c = (w>wmodel['blue'])&(w<wmodel['red'])
        mod = interpolate.splev(w[c],model)
        mod /= mod.mean()
        return (d[c]-mod)/abs(sky[c]+mod)**0.5

    # For `by hand' fitting
    if wsolution is None:
        import id_spec
        wsolution = id_spec.id_spec(sky,wmodel)
        x0 = numpy.arange(sky.size)
        w = sf.genfunc(x0,0.,wsolution)
        wsolution = sf.lsqfit(numpy.array([x,w]).T,'polynomial',3)
    else:
        import pylab
        tmpD = sky.copy()
        tmpD /= sky.mean()
        sky /= sky.mean()**2

        T = FFTFilter(tmpD,100)
        Dlines = spf.get_lines(x,T,nstd=7.)
        Dwlines = sf.genfunc(Dlines,0.,wsolution)
        w = sf.genfunc(x,0.,wsolution)
        c = (w>wmodel['blue'])&(w<wmodel['red'])
        mod = interpolate.splev(w[c],wmodel['model'])
        Slines = spf.get_lines(w[c],mod,nstd=7.)
        matches = []
        for j in range(Dwlines.size):
            diff = abs(Dwlines[j]-Slines)
            if diff.min()<5.*wmodel['scale']:
                matches.append([Dlines[j],Slines[diff.argmin()]])
        wsolution = sf.lsqfit(numpy.asarray(matches),'polynomial',3)

        start = wsolution['coeff'].flatten()
        coeff,ier = optimize.leastsq(refine,start,(x,tmpD,sky,wmodel['model']),
                    epsfcn=1e-5,maxfev=10000,ftol=1e-13,xtol=1e-13)
        wsolution = {'coeff':numpy.atleast_2d(coeff).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,wsolution)
        c = (w>wmodel['blue'])&(w<wmodel['red'])
        mod = interpolate.splev(w[c],wmodel['model'])
        mod /= mod.mean()
        pylab.plot(w[c],mod)
        pylab.plot(w[c],tmpD[c])
        pylab.show()
    w = sf.genfunc(x,0.,wsolution)
    rwsoln = sf.lsqfit(numpy.array([w,x]).T,'polynomial',3)
    return wsolution,rwsoln


def resampleData(data,owgrid,low,high,ysoln,xsoln,rwsoln):
    dy = (high-low)*2
    ygrid = iT.coords((dy,owgrid.size))[0]
    xgrid = sf.genfunc(owgrid,0.,rwsoln).repeat(dy)
    xgrid = xgrid.reshape((owgrid.size,dy)).T
    xgrid = sf.genfunc(xgrid.flatten(),ygrid.flatten(),xsoln['forw'])
    xgrid = xgrid.reshape(ygrid.shape)
    c = numpy.array([ygrid,xgrid])
    ygrid = ndimage.map_coordinates(ysoln[2],c)
    c = numpy.array([ygrid+low*2,xgrid])
    X = ndimage.map_coordinates(rotcoords2[1],c)
    Y = ndimage.map_coordinates(rotcoords2[0],c)
    c = numpy.array([Y,X])
    d = data.repeat(2,0).repeat(2,1)
    img = ndimage.map_coordinates(d,c,order=5)
    return iT.resamp(img,2)

def coadd2dspec(files, output, ymin=21, ymax=260, xmin=0, xmax=-1, shift=-1):
    """
    Added by Kuang-Han Huang (2014/09/12)
    Co-add 2-D spectra given by the list of file names.
    Also handle their RMS component (added by quardrature).
    Could also apply a shift (in arcsec); in this case, there must be two
    input images and the 2nd image is shifted DOWNWARDS relative to the 1st 
    image.
    """
    h = pyfits.open(files[0], ignore_missing_end=True)
    hdr = h[0].header  # use the header from the first file
    if xmax < 0:
        xmax = h[0].data.shape[1]
    data = h[0].data[ymin-1:ymax,xmin:xmax]
    var = h[1].data[ymin-1:ymax,xmin:xmax]
    h.close()
    if shift > 0:
        assert len(files) == 2, "The shift feature only works for two input images."
        print "Apply a shift of %.3f arcsec to the 2nd image (%s) relative" % (shift,files[1])
        print "to the 1st image (%s); note that the shift is DOWNWARD!!" % files[0]
        shiftpix = int(numpy.round(shift / 0.144))  # number of pixels to shift
        data1 = h[0].data[ymin-1:ymax,xmin:xmax]
        var1 = h[1].data[ymin-1:ymax,xmin:xmax]
        ny, nx = data1.shape
        data = numpy.zeros((ny + shiftpix, nx))
        var = numpy.zeros((ny + shiftpix, nx))
        data[shiftpix:,xmin:xmax] = data1
        var[shiftpix:,xmin:xmax] = var1
        h2 = pyfits.open(files[1], ignore_missing_end=True)
        data[:-shiftpix,xmin:xmax] = data[:-shiftpix,xmin:xmax] + h2[0].data[ymin-1:ymax,xmin:xmax]
        var[:-shiftpix,xmin:xmax] = var[:-shiftpix,xmin:xmax] + h2[1].data[ymin-1:ymax,xmin:xmax]
    else:    
        for f in files[1:]:
            h2 = pyfits.open(f, ignore_missing_end=True)
            data = data + h2[0].data[ymin-1:ymax,xmin:xmax]
            var = var + h2[1].data[ymin-1:ymax,xmin:xmax]
            h2.close()
    for i in range(len(files)):
        hdr['coadd%d' % (i+1)] = (files[i], "")
    hdu0 = pyfits.PrimaryHDU(data=data, header=hdr)
    hdu1 = pyfits.ImageHDU(data=var, header=hdr)
    hdulist = pyfits.HDUList([hdu0, hdu1])
    hdulist.writeto(output, clobber=True)


