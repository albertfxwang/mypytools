import numpy
from scipy import ndimage
import special_functions as sf
import indexTricks as iT

def crsub(data,sky,niter=4,nsig=5.,contrast=2.,sigfrac=0.3):
    for i in range(niter):
        map = cr_reject(data,sky,nsig,contrast,sigfrac)
        inmask = (1.-10000.*map)*data
        med5 = ndimage.median_filter(inmask,7)
        med5 *= map
        data = (1.-map)*data + med5
    return data

def cr_reject(sub,sky,nsig=5.,contrast=2.,sigfrac=0.3):
    gain = 5.
    rn = 25.
    sub = sub.copy()
    med5 = ndimage.median_filter(sub,5)
    med5 += sky

    n = 2
    blksub = sub.repeat(n,0).repeat(n,1)

    blksub = ndimage.laplace(blksub)*-0.5
    blksub[blksub<0] = 0.
    sub2 = iT.resamp(blksub,n,False)

    med5[med5<=0.] = 0.0001

    noise = (med5*gain + rn**2)**0.5/gain
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

    fmap = ndimage.maximum_filter(gmap,3)*sigmap
    fmap[fmap<(nsig*sigfrac)] = 0
    fmap[fmap>0] = 1
    map = fmap.copy()
    del gmap,fmap,stars,sigmap,med3,med7

    return map

