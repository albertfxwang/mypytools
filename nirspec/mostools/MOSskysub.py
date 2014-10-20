"""
Background subtract, resample, and coadd 2d spectra.
"""
import time
from mostools import spectools,correct_telluric,skysub
from special_functions import genfunc

import pyfits,scipy
from scipy import optimize,interpolate,ndimage,signal,stats


RESAMPLE = 0        # 1 to only resample
magicnum = -2**15



def offsety(data,outshape,offset):
    coords = spectools.array_coords(outshape)
    coords[0] += offset
    return ndimage.map_coordinates(data,coords,output=scipy.float64,cval=-999)


"""
Many of the pixels in the slit stack will be unused border pixels. This avoids
  taking the median of all of those (and therefore saves a bit of time).
"""
def fastmed(arr):
    tmp = arr.copy()
    tmp[scipy.isnan(tmp)] = 0.
    tmp = tmp.sum(0)
    tmp = tmp.sum(0)
    ends = scipy.where(tmp!=0)[0]
    left = ends.min()
    right = ends.max()+1
    tmp = stats.stats.nanmedian(arr[:,:,left:right],0)
    out = arr[0]*scipy.nan
    out[:,left:right] = tmp.copy()
    return out


def cr_reject(sub,sky,nsig=4.,contrast=1.):
    sub = sub.copy()
    med5 = ndimage.median_filter(sub,5)
    med5 += sky

    n = 2
    blksub = scipy.empty((sub.shape[0]*n,sub.shape[1]*n))
    yn,xn = sub.shape
    blksub = sub.repeat(n,0).repeat(n,1)

    blksub = ndimage.laplace(blksub)*-1.
    blksub[blksub<0] = 0.
    sumd = signal.correlate(blksub,scipy.ones((2,2)),mode='valid')/4.
    sub2 = sumd[::2,::2]
    sub2 /= 2.

    med5[med5<=0.] = 0.01

    noise = med5**0.5
    sigmap = sub2/noise
    med5 = ndimage.median_filter(sigmap,5)
    sigmap -= med5

    map = sigmap.copy()
    map[map<nsig] = 0
    map[map>0] = 1

    med3 = ndimage.median_filter(sub,3)
    med7 = ndimage.median_filter(med3,7)

    med3 -= med7
    med3 /= noise

    med3[med3<0.01] = 0.01

    stars = (map*sigmap)/med3

    stars[stars<contrast] = 0
    stars[stars>0.] = 1

    map *= stars

    gmap = ndimage.maximum_filter(map,3)*sigmap
    gmap[gmap<nsig] = 0
    gmap[gmap>0] = 1
    map = gmap.copy()

    return map




"""
doskysub()
"""
def doskysub(straight,ylen,xlen,sci,yback,sky2x,sky2y,ccd2wave,disp,mswave,offsets,cutoff,airmass):
    sci = sci.copy()

    # If cutoff is not a float, we are using the blueside
    if type(cutoff)==type([]):
        locutoff,hicutoff = cutoff
    else:
        locutoff = cutoff
        hicutoff = 10400.

    nsci = sci.shape[0]
    width = sci.shape[2]

    # Perform telluric correction
    coords = spectools.array_coords(sci[0].shape)
    x = coords[1].flatten()
    y = coords[0].flatten()

    for k in range(nsci):
        continue
        w = genfunc(x,y,ccd2wave[k])
        telluric = correct_telluric.correct(w,airmass[k],disp)
#       sci[k] *= telluric.reshape(sci[k].shape)
#    del coords,x,y,telluric

    # Create arrays for output images
#    outcoords = spectools.array_coords((ylen,xlen))
    outcoords = scipy.indices((ylen,xlen)).astype(scipy.float64)
    outcoords[1] *= disp
    outcoords[1] += mswave - disp*xlen/2.
    xout = outcoords[1].flatten()
    yout = outcoords[0].flatten()

    out = scipy.zeros((nsci,ylen,xlen))

    fudge = scipy.ceil(abs(offsets).max())
    bgimage = scipy.zeros((nsci,ylen+fudge,xlen))
    varimage = bgimage.copy()

    bgcoords = spectools.array_coords((ylen+fudge,xlen))
    bgcoords[1] *= disp
    bgcoords[1] += mswave - disp*xlen/2.

    #
    # Cosmic Ray Rejection and Background Subtraction
    #
    yfit = yback.flatten()
    ycond = (yfit>straight-0.4)&(yfit<straight+ylen-0.6)

    coords = spectools.array_coords(yback.shape)
    xvals = coords[1].flatten()
    yvals = coords[0].flatten()

    ap_y = scipy.zeros(0)
    aper = scipy.zeros(0)
    for k in range(nsci):
        xfit = genfunc(xvals,yfit-straight,ccd2wave[k])
        zfit = sci[k].flatten()

        x = xfit[ycond]
        y = yfit[ycond]
        z = zfit[ycond]

        # The plus/minus 20 provides a better solution for the edges
        wavecond = (x>locutoff-20.)&(x<hicutoff+20.)
        x = x[wavecond]
        y = y[wavecond]
        z = z[wavecond]

        # If only resampling...
        if RESAMPLE==1:
            coords = outcoords.copy()
            samp_x = genfunc(xout,yout,sky2x[k])
            samp_y = genfunc(xout,yout,sky2y[k])
            coords[0] = samp_y.reshape(coords[0].shape)
            coords[1] = samp_x.reshape(coords[1].shape)
            out[k] = scipy.ndimage.map_coordinates(sci[k],coords,output=scipy.float64,order=5,cval=-32768)

            out[k][xout.reshape(coords[1].shape)<locutoff] = scipy.nan
            out[k][xout.reshape(coords[1].shape)>hicutoff] = scipy.nan
            out[k][out[k]==-32768] = scipy.nan
            continue

        print "Determining sky for image %d"%(k+1)
        bgfit = skysub.skysub(x,y,z,disp)
        print "Subtracting sky"

        background = zfit.copy()
        a = time.time()
        for indx in range(background.size):
            x0 = xfit[indx]
            y0 = yfit[indx]
            if x0<locutoff-10 or x0>hicutoff+10:
                background[indx] = scipy.nan
            else:
                background[indx] = interpolate.bisplev(x0,y0,bgfit)
        sub = zfit-background
        sub[scipy.isnan(sub)] = 0.
        sky = sub*0.
        sky[ycond] = sub[ycond]
        sky = sky.reshape(sci[k].shape)
        sub = sky.copy()

        background[scipy.isnan(background)] = 0.

        # Note that 2d filtering may flag very sharp source traces!
        """
        sub = sub.reshape(sci[k].shape)
        sky = ndimage.median_filter(sky,5)
        diff = sub-sky
        model = scipy.sqrt(background.reshape(sci[k].shape)+sky)
        crmask = scipy.where(diff>4.*model,diff,0.)
        sub -= crmask
        sci[k] -= crmask
        """

        a = time.time()
        map = cr_reject(sub,background.reshape(sci[k].shape))

        inmask = (1.-10000.*map)*sub
        med5 = ndimage.median_filter(inmask,5)
        med5 *= map
        sub1 = (1.-map)*sub + med5
        crs = sub-sub1
        sub = sub1

        # Create straightened slit
        coords = outcoords.copy()
        samp_x = genfunc(xout,yout,sky2x[k])
        samp_y = genfunc(xout,yout,sky2y[k])
        coords[0] = samp_y.reshape(coords[0].shape)
        coords[1] = samp_x.reshape(coords[1].shape)
        out[k] = scipy.ndimage.map_coordinates(sci[k]-crs,coords,output=scipy.float64,order=5,cval=magicnum)
        vartmp = out[k].copy()
        out[k][xout.reshape(coords[1].shape)<locutoff] = scipy.nan
        out[k][xout.reshape(coords[1].shape)>hicutoff] = scipy.nan
        out[k][out[k]==magicnum] = scipy.nan

        # Output bgsub image
        coords = bgcoords.copy()
        bgy = bgcoords[0].flatten()+offsets[k]
        bgx = bgcoords[1].flatten()
        samp_x = genfunc(bgx,bgy,sky2x[k])
        samp_y = genfunc(bgx,bgy,sky2y[k])
        coords[0] = samp_y.reshape(coords[0].shape)
        coords[1] = samp_x.reshape(coords[1].shape)

#        varimage[k] = scipy.ndimage.map_coordinates(sci[k],coords,output=scipy.float64,order=5,cval=magicnum)
        crs = scipy.ndimage.map_coordinates(crs,coords,output=scipy.float64,order=1,cval=magicnum)
        crs[crs>0.3] = scipy.nan
        varimage[k] = crs+vartmp

        # Only include good data (ie positive variance, wavelength
        #   greater than dichroic cutoff)
        cond = (bgcoords[0]+offsets[k]<0.)|(bgcoords[0]+offsets[k]>ylen)
        cond = (varimage[k]<0)|cond
        cond = (bgcoords[1]<locutoff)|(bgcoords[1]>hicutoff)|cond
        varimage[k][cond] = scipy.nan

        bgimage[k] = scipy.ndimage.map_coordinates(sub,coords,output=scipy.float64,order=5,cval=magicnum)
        bgimage[k][cond] = scipy.nan
        bgimage[k][bgimage[k]==magicnum] = scipy.nan # Shouldn't be
                                 #   necessary...

    if RESAMPLE==1:
        return out,bgimage,varimage,outcoords[1,0]
    if bgimage.shape[0]>1:    
        bgimage = fastmed(bgimage)
        varimage = fastmed(varimage)/nsci
    elif bgimage.ndim==3:
        bgimage = bgimage[0].copy()
        varimage = varimage[0].copy()


    return out,bgimage,varimage,outcoords[1,0]
