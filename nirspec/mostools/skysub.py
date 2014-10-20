import scipy,special_functions
from scipy import ndimage,interpolate

WIDE = 110
AGGRESSIVE = False

def skysub(x,y,z,scale):
    """
    skysub(x,y,z,scale)

    Routine to determine the 2d background from data. (x,y) are the
      coordinates of the data, usually in the *corrected* frame.

    Inputs:
      x     - 1d array describing x-coordinate, usually wavelength
      y     - 1d array describing y-coordinate, usually corrected spatial
            position
      z     - data each position (x,y)
      scale - approximate output scale (for knot placement). It is not, in
            general, possible to calculate this from x because the
            input coordinates are not on a regular grid.

    Outputs:
      2d spline model of the background
    """
    cond = (scipy.isfinite(z))&(z>0.)
    x = x[cond]
    y = y[cond]
    z = z[cond]

    x0 = x.copy()
    y0 = y.copy()
    z0 = z.copy()

    height = int(y.max()-y.min())
    width = int(x.max()-x.min())
    npoints = x.size

    midpt = y.mean()

    """
    Very wide slits need special attention. Here we fit a first order
      correction to the slit and subtract it away before doing the high
      pixel rejection (the problem is if there is a small gradient across
      a wide slit, the top and bottom pixels may differ significantly,
      but these pixels may be close in *wavelength* and so locally (on
      the CCD) low pixels will be rejected in the smoothing
    """
    if height>WIDE:
        zbak = z.copy()
        args = y.argsort()
        revargs = args.argsort()
        ymodel = ndimage.percentile_filter(z[args],30.,size=height)[revargs]
        fit = special_functions.lsqfit(ymodel,'polynomial',1)

        if fit['coeff'][1]*float(ymodel.size)/fit['coeff'][0]<0.05:
            pass
        else:
            ymodel = special_functions.genfunc(scipy.arange(ymodel.size),0,fit)
            ymodel -= ymodel.mean()
            z -= ymodel

    # Filter locally (in wavelength space) high points
    args = x.argsort()
    revargs = args.argsort()

    smooth = ndimage.percentile_filter(z[args],35.,size=height)[revargs]

    diff = z-smooth
    # We assume poisson statistics....
    var = scipy.sqrt(scipy.fabs(z))
    sigma = diff/var

    args = y.argsort()
    revargs = args.argsort()

    t = ndimage.median_filter(sigma[args],9)
    t = ndimage.gaussian_filter(t,width)#[revargs]
    # Source detection/rejection
    # Reject yvalues > 1. sigma, and weight remaining pixels
    w = (1.0-t)/abs(z[args])

    if AGGRESSIVE:
        g = scipy.where(w<=0,0,1)
        g = ndimage.maximum_filter(g,width*3)
        g = ndimage.minimum_filter(g,width*7)

        s = sigma[args].copy()
        b = ndimage.minimum_filter(g,width*5)
        xi = scipy.arange(t.size)
        fitdata = scipy.empty((xi[g==1].size,2))
        fitdata[:,0] = xi[g==1].copy()
        fitdata[:,1] = t[g==1].copy()
        fit = special_functions.lsqfit(fitdata,'polynomial',3)
        fit = special_functions.genfunc(xi,0.,fit)

        diff = (t-fit)[b==1]
        s = diff.std()
        while (abs(t-fit)[(g==1)&(b==0)]>2.5*s).any():
            g = b.copy()
            b = ndimage.minimum_filter(g,width*5)
            fitdata = scipy.empty((xi[g==1].size,2))
            fitdata[:,0] = xi[g==1].copy()
            fitdata[:,1] = t[g==1].copy()
            fit = special_functions.lsqfit(fitdata,'polynomial',3)
            fit = special_functions.genfunc(xi,0.,fit)

            diff = (t-fit)[b==1]
            s = diff.std()

        w *= g

    skycond = ((w>0.)&(z>0))
    x = x[skycond]
    y = y[skycond]
    z = z[skycond]

    # Reject residual high pixels (and very low pixels too!)
    args = x.argsort()
    revargs = args.argsort()
    smooth = ndimage.median_filter(z[args],height/4.)[revargs]

    diff = z-smooth
    var = scipy.sqrt(smooth)

    cond = abs(diff)<4.*var
    x = x[cond]
    y = y[cond]
    z = z[cond]

    kx = 3
    ky = 1


    # If the slit is long, return to original data and increase the order
    #   of the y-fit.
    if height>WIDE:
        z = zbak[skycond]
        z = z[cond].astype(scipy.float64)

    #    if height>WIDE*1.5:
    #        ky = 3

        cond = z>0.
        x = x[cond]
        y = y[cond]
        z = z[cond]

    w = 1./z

    if x.size<5.*width:
        kx = 1
        ky = 1

    # Create knots...
    
    innertx = scipy.arange(x.min()+scale/2.,x.max()-scale/2.,3.*scale/4.)
    """
    tx = scipy.zeros(innertx.size+kx*2+2)
    tx[0:kx+1] = x.min()
    tx[kx+1:innertx.size+kx+1] = innertx.copy()
    tx[innertx.size+kx+1:] = x.max()
    """
    tx = scipy.linspace(x.min(),x.max(),innertx.size)
    xsort = scipy.sort(x)
    tmp = [x.min()]
    num = []
    cnt = 0
    j = 1
    for i in range(xsort.size):
        while xsort[i]>tx[j]:
            if cnt>0:
                if len(num)==0 or cnt>1 or num[-1]>1:
                    tmp.append(tx[j])
                    num.append(cnt)
                    cnt = 0
            j += 1
        cnt += 1
    tmp.append(x.max())
    tx = scipy.asarray(tmp)
    ty = scipy.zeros(ky*2+2)
    ty[0:ky+1] = y.min()
    ty[ky+1:] = y.max()

    #del innertx
    # ...and fit.
    bgfit = interpolate.bisplrep(x,y,z,w,tx=tx,ty=ty,kx=kx,ky=ky,task=-1,nxest=tx.size,nyest=ty.size)
    del x,y,z,w,tx,ty
    return bgfit

