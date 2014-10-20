import scipy
import special_functions as sf


def new_xsoln(data,xord=3,yord=3,tol=2.):
    data = data.copy()

    height,width = data.shape

    indx = height/4
    mid = data[indx]

    xvals = scipy.arange(mid.size)
    middle_lines = get_lines(xvals,mid,False)

    ref = []
    for row in range(height):
        ref.append(get_lines(xvals,data[row],False))
    lines = []
    for l in middle_lines:
        """
        off is used to track the systematic shift of lines as we move
            up and down the slit. This is necessary for strongly
            curved lines that might be substantially shifted wrt
            the central row.
        """
        off = 0.
        for row in range(indx,height):    # The top of the slit
            if len(ref[row])==0:
                continue
            p = ref[row]+off
            diff = abs(l-p)
            if diff.min()<tol:
                pos = diff.argmin()
                lines.append([ref[row][pos],row,l])
                off = l-ref[row][pos]
        off = 0.
        for row in range(indx-1,-1,-1):    # The bottom of the slit
            if len(ref[row])==0:
                continue
            p = ref[row]+off
            diff = abs(l-p)
            if diff.min()<tol:
                pos = diff.argmin()
                lines.append([ref[row][pos],row,l])
                off = l-ref[row][pos]
    lines = scipy.asarray(lines)
    soln = sf.lsqfit(lines,'chebyshev',xord,yord)

    tmp = lines[:,0].copy()
    lines[:,0] = lines[:,2].copy()
    lines[:,2] = tmp.copy()
    soln2 = sf.lsqfit(lines,'chebyshev',xord,yord)

    return {'back':soln,'forw':soln2}


def update_xsoln(data,soln,xord=3,yord=3,skip=False):
    data = data.copy()

    height,width = data.shape

    if skip:
        slice = data.mean(1)
        indx = slice.argsort()[slice.size/4]
        mid = data[indx]
        badrow = scipy.where(slice>slice[indx]+5.*(slice[indx])**0.5)[0]
        if badrow.size==0:
            badrow = scipy.array([-1])
    else:
        indx = height/2
        mid = data[indx]
        badrow = scipy.array([-1])

    xvals = scipy.arange(mid.size)
    middle_lines = get_lines(xvals,mid,False)
    straight = sf.genfunc(middle_lines,indx,soln['back'])
    ref = []
    lines = []

    for row in range(height):
        if abs(row-badrow).min()==0:
            continue
        current =  get_lines(xvals,data[row],False)
        if current.size==0:
            continue
        guess = sf.genfunc(current,row,soln['back'])
        for l in straight:
            diff = abs(l-guess)
            if diff.min()<2:
                pos = diff.argmin()
                lines.append([current[pos],row,l])
    lines = scipy.asarray(lines)
    newsoln = sf.lsqfit(lines,'chebyshev',xord,yord)

    tmp = lines[:,0].copy()
    lines[:,0] = lines[:,2].copy()
    lines[:,2] = tmp.copy()
    del tmp
    newsoln2 = sf.lsqfit(lines,'chebyshev',xord,yord)

    return {'back':newsoln,'forw':newsoln2}


def newModMatch(data,model,scale,order,center,sparse=False):
    import numpy
    import special_functions as sf
    from scipy import ndimage,interpolate
    x0 = numpy.arange(data.size).astype(numpy.float32)
    datalines = get_lines(x0,data)
    Xstart = int(datalines.min()-15.)
    Xend = int(datalines.max()+15.)
    if Xstart<0:
        Xstart = 0
    if Xend>data.size:
        Xend = data.size
    tdata = data[Xstart:Xend].copy()
    wave = numpy.linspace(center-scale*data.size/2,center+scale*data.size/2,data.size)[Xstart:Xend].copy()
    cond = wave<10300.
    mod = interpolate.splev(wave,model['matched'])
    fftd = numpy.fft.rfft(ndimage.gaussian_filter(tdata[cond],17))
    fftm = numpy.fft.rfft(ndimage.gaussian_filter(mod[cond],17)[::-1])
    T = numpy.fft.fftshift(numpy.fft.irfft(fftd*fftm)).real
    offset = T.size/2-T.argmax()
    Wstart = wave[offset]

    pars = numpy.zeros(order+1)
    pars[0] = Wstart-Xstart*scale
    pars[1] = scale
    cov = numpy.zeros(order+1)
    cov[0] = 1.
    cov[1] = 0.05
    for i in range(2,order+1):
        cov[i] = (1./data.size)**i

    coeff = pars
    def optfunc(pars,x,d,mod,sig):
        if numpy.isnan(pars).any():
            return -1e200
        coeff = numpy.atleast_2d(numpy.array(pars)).T
        fit = {'coeff':coeff,'type':'polynomial'}
        x0 = sf.genfunc(x,0.,fit)
        try:
            model = interpolate.splev(x0,mod)
        except:
            return -1e200
        model /= model.mean()
        resid = (model-d)/sig
        return -0.5*(resid**2).sum()

    spec = ndimage.gaussian_filter(tdata,7)
    spec /= spec.mean()
    sig = spec**0.5
    logp = optfunc(pars,x0[Xstart:Xend],spec,model['wide'],sig)
    niter = 2000
    for indx in xrange(niter):
        prop = coeff+numpy.random.randn(cov.size)*cov*10**(numpy.random.random(cov.size)*2.3-2.)
        tlogp = optfunc(prop,x0[Xstart:Xend],spec,model['wide'],sig)
        accept = 0
        if tlogp>logp:
            accept = 1
        else:
            if tlogp-logp>numpy.log(numpy.random.random()):
                accept = 1
        if accept==1:
            logp = tlogp
            coeff = prop.copy()

    spec = tdata.copy()
    spec /= spec.mean()
    sig = spec**0.5
    logp = optfunc(pars,x0[Xstart:Xend],spec,model['matched'],sig)
    c = numpy.empty((niter,coeff.size))
    for indx in xrange(niter):
        prop = coeff+numpy.random.randn(cov.size)*cov*10**(numpy.random.random(cov.size)*2.3-2.)
        tlogp = optfunc(prop,x0[Xstart:Xend],spec,model['matched'],sig)
        accept = 0 
        if tlogp>logp:
            accept = 1 
        else:
            if tlogp-logp>numpy.log(numpy.random.random()):
                accept = 1 
        if accept==1:
            logp = tlogp
            coeff = prop.copy()
        c[indx] = coeff.copy()
    for indx in range(cov.size):
        cov[indx] = c[:,indx].std()
    mlogp = logp
    mcoeff = coeff.copy()
    for indx in xrange(niter):
        prop = coeff+numpy.random.randn(cov.size)*cov*10**(numpy.random.random(cov.size)*2.3-2.)
        tlogp = optfunc(prop,x0[Xstart:Xend],spec,model['matched'],sig)
        accept = 0
        if tlogp>logp:
            accept = 1
        else:
            if tlogp-logp>numpy.log(numpy.random.random()):
                accept = 1
        if accept==1:
            logp = tlogp
            coeff = prop.copy()
            if logp>mlogp:
                mlogp = logp
                mcoeff = coeff.copy()
    fit = {'coeff':scipy.atleast_2d(mcoeff).T,'type':'polynomial'}
    w = sf.genfunc(x0[Xstart:Xend],0.,fit)
    m = interpolate.splev(w,model['matched'])
    cond = (w<10400)
    m = m[cond]
    m /= m.mean()
    import pylab
    pylab.plot(w[cond],m)
    pylab.plot(w[cond],spec[cond])
    pylab.show()
    return fit




#def newModMatch(data,model,scale,order,center,sparse=False):
    import numpy,pylab
    from scipy import interpolate,optimize,ndimage

    x0 = numpy.arange(data.size)

    datalines = get_lines(x0,data)
    if datalines.size<15:
        sparse = True
    Xstart = int(datalines.min()-15.)
    Xend = int(datalines.max()+15.)
    if Xstart<0:
        Xstart = 0
    if Xend>data.size:
        Xend = data.size
    tdata = data[Xstart:Xend].copy()
    wave = numpy.linspace(center-scale*tdata.size/2,center+scale*tdata.size/2,tdata.size)
    cond = (wave<10300)
    mod = interpolate.splev(wave,model['matched'])
    fftd = numpy.fft.rfft(ndimage.gaussian_filter(tdata[cond],17))
    fftm = numpy.fft.rfft(ndimage.gaussian_filter(mod[cond],17)[::-1])
    fftd[100:] = 0.
    fftm[100:] = 0.
    T = numpy.fft.fftshift(numpy.fft.irfft(fftd*fftm)).real
    offset = T.size/2-T.argmax()
#    offset = T.size/2 - 1747
    Wstart = wave[offset]
    print "Assuming the starting wavelength is %6.1f"%(Wstart)
#    pylab.plot(T)
#    pylab.show()
    if sparse is True:
        nlines = datalines.size
        n = order+2
        c = numpy.atleast_2d(numpy.array([Wstart,scale])).T
        fit = {'coeff':c,'type':'polynomial'}
        first = True
        while n<nlines:
            end = datalines[n]+15.
            if end>data.size:
                end = data.size
            tdata = data[Xstart:end].copy()
            x = x0[Xstart:end].copy()
            if first is True:
                W = sf.genfunc(x-Xstart,0.,fit)
            else:
                W = sf.genfunc(x,0.,fit)
            cond = (W<10400)
            W = W[cond]
            tdata = tdata[cond]
            x = x[cond]
            mod = interpolate.splev(W,model['matched'])
            fftd = numpy.fft.rfft(tdata)
            fftm = numpy.fft.rfft(mod[::-1])
            T = numpy.fft.fftshift(numpy.fft.irfft(fftd*fftm)).real
            offset = T.size/2-T.argmax()
            Wlines = get_lines(x-offset,mod)
            Wwave = get_lines(W,mod)
            matches = []
            for i in range(Wlines.size):
                diff = abs(Wlines[i]-datalines[:n])
                if diff.min()<2:
                    matches.append([datalines[diff.argmin()],Wwave[i]])
            if len(matches)<order+1:
                n += 1
                continue
            matches = numpy.asarray(matches)
            fit = sf.lsqfit(matches,'polynomial',order)
            n += 1
            first = False

        return fit

    coeff = numpy.ones(order+1)
    rescale = numpy.array([Wstart-Xstart*scale,scale])
    for i in range(order-1):
        rescale = numpy.append(rescale,1./(data.size/2)**(i+2))

    def optFunc(p,x,spec,mod):
        if (numpy.isnan(p)).any() or p[0]<0.5 or p[1]<0.8 or p[1]>1.2:
            return spec
        coeff = numpy.atleast_2d(numpy.array(p)*rescale).T
        fit = {'coeff':coeff,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,mod)
        cond = (w<9800.)
        m = m[cond]
        m /= m.mean()
        resid = (spec[cond]-m)/abs(spec[cond]+m)**0.5
        return resid/resid.size

    import pylab
    from scipy import ndimage
    n = 2*order+2
    end = datalines[n]+15
    if end<Xstart+(Xend-Xstart)/4:
        end = Xstart + (Xend-Xstart)/4
    widespec = ndimage.gaussian_filter(data,15)
    x = x0[Xstart:end].copy()
    spec = widespec[Xstart:end].copy()
    spec /= spec.mean()
    coeff,ier = optimize.leastsq(optFunc,coeff,(x,spec,model['wide']),factor=99.,ftol=1e-12,xtol=1e-12,epsfcn=1e-11)
    while end<Xend:
        x = x0[Xstart:end].copy()
#        spec = widespec[Xstart:end].copy()
#        spec /= spec.mean()
#        coeff,ier = optimize.leastsq(optFunc,coeff,(x,spec,model['wide']),factor=99.,ftol=1e-12,xtol=1e-12,epsfcn=1e-11)
        spec = data[Xstart:end].copy()
        spec /= spec.mean()
        coeff,ier = optimize.leastsq(optFunc,coeff,(x,spec,model['matched']),factor=99.,ftol=1e-12,xtol=1e-12,epsfcn=1e-11)

        fit = {'coeff':scipy.atleast_2d(coeff*rescale).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model['matched'])
        cond = (w<10400)
        m = m[cond]
        m /= m.mean()
        coeff,ier = optimize.leastsq(optFunc,coeff,(x,spec,model['matched']),factor=99.,ftol=1e-12,xtol=1e-12,epsfcn=1e-11)
        end += (Xend-Xstart)/10.
    x = x0[Xstart:Xend].copy()
    spec = data[Xstart:Xend].copy()
    spec /= spec.mean()

    coeff,ier = optimize.leastsq(optFunc,coeff,(x,spec,model['matched']),
                    ftol=1e-12,xtol=1e-12,epsfcn=1e-7)
    fit = {'coeff':scipy.atleast_2d(coeff*rescale).T,'type':'polynomial'}
    w = sf.genfunc(x,0.,fit)
    m = interpolate.splev(w,model['matched'])
    cond = (w<10400)
    m = m[cond]
    m /= m.mean()
    pylab.plot(w[cond],m)
    pylab.plot(w[cond],spec[cond])
    pylab.show()
    return fit

    w = sf.genfunc(x,0.,fit)
    m = interpolate.splev(w,model['matched'])
    rhs = m/err
    fit = numpy.array(numpy.linalg.lstsq(op,rhs)[0])
    d = numpy.dot(linmod.T,fit)
    pylab.plot(w,m)
    pylab.plot(w,d)

    w0 = w.copy()
    fit = {'coeff':scipy.atleast_2d(pars).T,'type':'polynomial'}
    w = sf.genfunc(x,0.,fit)
    m = interpolate.splev(w,model['matched'])
    rhs = m/err
    fit = numpy.array(numpy.linalg.lstsq(op,rhs)[0])
    d2 = numpy.dot(linmod.T,fit)
    pylab.plot(w,d2)

    print (w0-w).std()
    pylab.figure()
    pylab.plot(w0-w)
    pylab.show()

    return fit


    df


def modelmatch(data,model,scale,order,fitrange,minwave=None,maxwave=None):
    """
    Match a spectrum against a model of the spectrum.
    """
    from scipy import ndimage,stats,interpolate,signal,optimize
    data = data.copy()

    if data.ndim==2:
        spec = stats.stats.nanmedian(data,0)
    else:
        spec = data.copy()

    spec[scipy.isnan(spec)] = 0.
    spec[spec<0.] = 0.
    spec[spec==0] = scipy.median(spec)
    wide = ndimage.gaussian_filter(spec,5.)

    lowave,hiwave = fitrange
    if lowave is None:
        lowave = 5150.
    if hiwave is None:
        hiwave = 9000.
    if hiwave>10400.:
        hiwave = 10400.
    if maxwave is None:
        maxwave = 10400.
    if maxwave>10400.:
        maxwave = 10400.
    if hiwave>maxwave:
        hiwave = maxwave
    if minwave is None:
        minwave = 3000.
    if lowave<minwave:
        lowave = minwave

    outw = scipy.arange(minwave,hiwave,scale)
    out = outw*0.
    outc = out*0.

    x = scipy.arange(wide.size)
    l = get_lines(x,spec,nstd=10.)
    l = l[wide[l.astype(scipy.int16)]>scipy.median(wide)]
    offset = int(l.min()-50.)

    if offset<0:
        offset = 0.
    xspec = wide[offset:]

    """
    `width' is the range over which the wavelength solution is assumed to be
        linear.
    """
    width = 250.
    w0 = lowave
    while w0<hiwave:
        w = scipy.arange(w0,w0+width,scale)
        mod = interpolate.splev(w,model['wide'])
        corr = signal.correlate(mod,xspec,'valid')

        for j in range(corr.size):
            x = xspec[j:j+mod.size]
            med = scipy.median(x/mod)
            corr[j] /= ((mod*med-x)**2/x).sum()

        corr /= corr.mean()
        corr = ndimage.gaussian_filter(corr,5)
        w = w0-scipy.arange(corr.size)*scale
        w = w[::-1]
        corr = corr[::-1]
        cond = (outw>=w[0])&(outw<=w[-1])
        mod = interpolate.splrep(w,corr,s=0)
        out[cond] += interpolate.splev(outw[cond],mod)
        outc[cond] += 1
        w0 += width

    out[outc>0] /= outc[outc>0]
    start = outw[out.argmax()]-offset*scale

    def dofit(p,x,data,model):
        if scipy.isnan(p).any():
            return x*0.+1e7
        if abs(p[2]-scale)/scale>0.15:
            return x*0.+1e7
        fit = {'coeff':scipy.atleast_2d(p[1:]).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model)
        m *= p[0]
        chi = (m-data)/abs(data)**0.5
        cond = ~scipy.isnan(chi)
        cond = cond&scipy.isfinite(chi)
        cond = cond&(w>lowave)&(w<maxwave)
        badval = abs(chi[cond]).max()
        chi[~cond] = 2*badval
        return chi
#        if chi[cond].size<x.size/2:
#            return x*0.+1e7
#        print chi[cond].size,w
#        return chi[cond]/chi[cond].size

    spec = spec[offset:]
    wide = wide[offset:]

    w = scipy.arange(start+offset*scale,11400.,scale)
    m = interpolate.splev(w,model['matched'])[:spec.size]
    m[w[:spec.size]>10400.] = 0.
    print offset,start
    x = scipy.arange(spec.size).astype(scipy.float32)+offset

#    import pylab
#    pylab.plot(spec)
#    pylab.show()
    print spec.size,m.size
    initial_pars = [scipy.median(spec/m),start,scale]
    diag = [1.,1./start,1./scale]
    for i in range(order+2-len(initial_pars)):
        initial_pars.append(0.)
        diag.append(1e4**(i+2))

#    import pylab
#    pylab.plot(spec)
#    pylab.plot(m*initial_pars[0])
#    pylab.show()

    # Pixels with very low values are assumed invaled
    cond = spec>1e-5
    spec = spec[cond]
    wide = wide[cond]
    x = x[cond]

    import numpy
#    diag = numpy.ones(len(diag))
    coeff,ier = optimize.leastsq(dofit,initial_pars,(x,wide,model['wide']),
                        maxfev=1e5,epsfcn=1e-15,diag=diag,ftol=1e-15,factor=99.)
    coeff,ier = optimize.leastsq(dofit,coeff,(x,spec,model['matched']),
                        maxfev=1e5,epsfcn=1e-15,diag=diag,ftol=1e-15,factor=99.)
    coeff,ier = optimize.leastsq(dofit,coeff,(x,spec,model['matched']),
                        maxfev=1e5,epsfcn=1e-5,diag=diag,ftol=1e-15,factor=99.)

    fit = {'coeff':scipy.atleast_2d(coeff[1:]).T,'type':'polynomial'}
    if 1==1: #DEBUG
        import pylab
        print coeff
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model['matched'])*coeff[0]
        pylab.plot(w,m)
        pylab.plot(w,spec)
        pylab.show()

    return fit


def skywave(data,skymodel,scale,order=3,cutoff=None,hicutoff=9100.,minwave=None):
    from scipy import ndimage,stats,interpolate,signal,optimize
    import pylab

    data = data.copy()

    if data.ndim==2:
        sky = stats.stats.nanmedian(data,0)
    else:
        sky = data.copy()
    sky[scipy.isnan(sky)] = 0.
    sky[sky<0.] = 0.
    sky[sky==0.] = scipy.median(sky)
    osky = sky.copy()
    sky = ndimage.gaussian_filter(sky,5.)

    if cutoff is None:
        outw = scipy.arange(3000.,9000.,1.)
    else:
        outw = scipy.arange(cutoff-0.5*osky.size*scale,9000.,1.)
    if minwave is not None:
        outw = outw[outw>minwave]

    out = outw*0.
    outc = out*0.

    x = scipy.arange(osky.size)
    l = get_lines(x,sky,nstd=15.)

    l = l[sky[l.astype(scipy.int16)]>scipy.median(sky)]
    offset = int(l.min()-50.)
    if offset<0:
        offset = 0
    sky = sky[offset:]

    width = 250.
    if cutoff is None:
        w0 = 5150.-width
    else:
        w0 = cutoff-width

    while w0+width<hicutoff:
        w0 += width
        if cutoff is not None and w0<cutoff:
            continue
        w1 = scipy.arange(w0,w0+width,scale)
        t1 = interpolate.splev(w1,skymodel['wide'])
        corr = signal.correlate(t1,sky,'valid')

        chi = corr*0.
        for j in range(corr.size):
            med = scipy.median(sky[j:j+t1.size]/t1)
            chi[j] = ((t1*med-sky[j:j+t1.size])**2/sky[j:j+t1.size]).sum()

        ratio = corr.copy()/chi.copy()
        ratio /= ratio.mean()
        ratio = ndimage.gaussian_filter(ratio,5)
        mywave = w0-scipy.arange(ratio.size)*scale

        mywave = mywave[::-1]
        ratio = ratio[::-1]
        cond = (outw>=mywave[0])&(outw<=mywave[-1])

        mod = interpolate.splrep(mywave,ratio,s=0)
        out[cond] += interpolate.splev(outw[cond],mod)
        outc[cond] += 1.

    out[outc>0] /= outc[outc>0]
    start = outw[out.argmax()]-offset*scale
    sky = osky[offset:].copy()
    sky_wide = ndimage.gaussian_filter(sky,5)


    cutoff = 5500.
    """ Optimization function for refining the wavelength solution. """
    def dofit(p,x,data,model):
        if scipy.isnan(p).any():
            return x*0.+1e7
        fit = {'coeff':scipy.atleast_2d(p[1:]).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model)
        m *= p[0]
        chi = (m-data)/abs(data)**0.5
        cond = ~scipy.isnan(chi)
        cond = cond&scipy.isfinite(chi)
        cond = cond&(w>cutoff)&(w<10400.)
        return chi[cond]/chi[cond].size

    w = scipy.arange(start+offset*scale,10400.,scale)
    m = interpolate.splev(w,skymodel['matched'])

    x = scipy.arange(sky.size).astype(scipy.float32)+offset
    initial_pars = [scipy.median(sky/m[:sky.size]),start,scale]
    diag = [1.,1e3,1.]
    for i in range(order+2-len(initial_pars)):
#        initial_pars.append(1e-4**(i+2))
        initial_pars.append(0.)
        diag.append(1e3**(i+2))

    cond = sky>1e-5
    sky = sky[cond]
    sky_wide = sky_wide[cond]
    x = x[cond]

    coeff,ier = optimize.leastsq(dofit,initial_pars,
                        (x,sky_wide,skymodel['wide']),maxfev=1e5,epsfcn=1e-15,diag=diag,xtol=1e-15,ftol=1e-15)
    coeff,ier = optimize.leastsq(dofit,coeff,
                        (x,sky,skymodel['matched']),maxfev=1e5,epsfcn=1e-15,diag=diag,xtol=1e-15,ftol=1e-15)


    return {'coeff':scipy.atleast_2d(coeff[1:]).T,'type':'polynomial'}


def arcwave(sky,arc,arcmodel,skymodel,scale,order):
    from scipy import ndimage,stats,interpolate,optimize
    sky = sky.copy()
    arc = arc.copy()

    sky = scipy.median(sky,0)

    x = scipy.arange(sky.size)
    x_orig = x.copy()

    wave = scipy.arange(3000.,10000.,scale)
    arc_wide = ndimage.gaussian_filter(arc,5)
    m = interpolate.splev(wave,arcmodel['norm'])

    a = arc.copy()
    aw = arc_wide.copy()
    arc = a[:a.size/2.]

    x = x_orig[:a.size/2.]
    arclines = get_lines(x,arc)
    fit = scipy.zeros(3*arclines.size+1)
    index = 1
    for i in range(arclines.size):
        fit[index] = 1.
        fit[index+1] = arclines[i]
        fit[index+2] = 15.*scale
        index += 3
    arc_wide = sf.ngauss(x,fit)
    """
    Do an approximate chi-square between the sky model and the data over a
        range of offsets using a broadened data and sky model.
    """
    max = 0.
    mid = 0

    delta = scale/10.
    s = scipy.arange(scale-delta,scale+delta,delta/10.)
    for stmp in s:
        wtmp = scipy.arange(2000.,10000.,stmp)
        m = interpolate.splev(wtmp,arcmodel['norm'])
        conv = scipy.empty(m.size-arc_wide.size+1)
        for i in range(conv.size):
            tmp = m[i:i+arc_wide.size].copy()
            if tmp.max()<0.1:
                conv[i] = 0.
                continue
            conv[i] = (tmp*arc_wide).sum()
            conv[i] = 1./((tmp-arc_wide)**2).sum()
        curr = conv.max()
        if curr>max:
            mid = conv.argmax()
            scale = stmp
            max = conv.max()
            wave = wtmp.copy()

    """
    Refine the starting wavelength position using the 'true' (ie narrow) model
        of the sky. Searches for a minimum around the minimum found in the
        previous optimization.
    """
    m = interpolate.splev(wave,arcmodel['matched'])
    conv = scipy.empty(m.size-arc.size+1)
    for i in range(conv.size):
        tmp = m[i:i+arc.size].copy()
        ratio = arc.max()/tmp.max()
        if tmp.max()<1.:
            conv[i] = 0.
            continue
        tmp *= ratio
        conv[i] = (tmp*arc).sum()
    pt = conv[mid-50:mid+51].argmax()+mid-50


    initial_pars = [wave[pt],scale]
    for i in range(order+1-len(initial_pars)):
        initial_pars.append(0.)
    modellines = get_lines(wave,m,std=10.)
    modellines = modellines[modellines>wave[pt]]
    modellines = arcmodel['lines']
    modellines = modellines[modellines>wave[pt]]

    fit = {'coeff':scipy.atleast_2d(initial_pars).T,'type':'polynomial'}

    for o in [1,2]:
        w = sf.genfunc(arclines,0.,fit)
        matches = []
        for j in range(w.size):
            diff = abs(w[j]-modellines)
            if diff.min()<5.*scale:
                matches.append([arclines[j],modellines[diff.argmin()]])
        fit = sf.lsqfit(scipy.asarray(matches),'polynomial',o)

    left_matches = [i for i in matches]    
    wmin = sf.genfunc(a.size*0.45,0.,fit)

    arc = a[a.size/2.:].copy()
    x = scipy.arange(arc.size).astype(scipy.float32)+a.size/2.
    arclines = get_lines(x,arc)
    fit = scipy.zeros(3*arclines.size+1)
    index = 1
    for i in range(arclines.size):
        fit[index] = 1.
        fit[index+1] = arclines[i]
        fit[index+2] = 10.*scale
        index += 3
    arc_wide = sf.ngauss(x,fit)
    """
    Do an approximate chi-square between the sky model and the data over a
        range of offsets using a broadened data and sky model.
    """
    max = 0.
    mid = 0
    delta = scale/10.
    s = scipy.arange(scale-delta,scale+delta,delta/10.)
    for stmp in s:
        wtmp = scipy.arange(wmin,10000.,stmp)
        m = interpolate.splev(wtmp,arcmodel['norm'])
        conv = scipy.empty(m.size-arc_wide.size+1)
        for i in range(conv.size):
            tmp = m[i:i+arc_wide.size].copy()
            if tmp.max()<0.1:
                conv[i] = 0.
                continue
            conv[i] = (tmp*arc_wide).sum()
        curr = conv.max()
        if curr>max:
            mid = conv.argmax()
            scale = stmp
            max = conv.max()
            wave = wtmp.copy()
    """
    Refine the starting wavelength position using the 'true' (ie narrow) model
        of the sky. Searches for a minimum around the minimum found in the
        previous optimization.
    """
    m = interpolate.splev(wave,arcmodel['matched'])
    conv = scipy.empty(m.size-arc.size+1)
    for i in range(conv.size):
        tmp = m[i:i+arc.size].copy()
        ratio = arc.max()/tmp.max()
        if tmp.max()<1.:
            conv[i] = 0.
            continue
        tmp *= ratio
        conv[i] = (tmp*arc).sum()
    pt = conv[mid-50:mid+51].argmax()+mid-50
    wavept = wave[pt]

    initial_pars = [wavept,scale]
    for i in range(order+1-len(initial_pars)):
        initial_pars.append(0.)
    modellines = get_lines(wave,m,std=10.)
    modellines = modellines[modellines>wavept]
    modellines = arcmodel['lines']
    modellines = modellines[modellines>wavept]

    fit = {'coeff':scipy.atleast_2d(initial_pars).T,'type':'polynomial'}
    for o in [1,2]:
        # The (o-2) bit is to correct the offset after the first loop
        w = sf.genfunc(arclines+(o-2)*a.size/2.,0.,fit)
        matches = []
        for j in range(w.size):
            diff = abs(w[j]-modellines)
            if diff.min()<5.*scale:
                matches.append([arclines[j],modellines[diff.argmin()]])
        fit = sf.lsqfit(scipy.asarray(matches),'polynomial',o)

    arc = a.copy()
    arc_wide = aw.copy()

    w = sf.genfunc(arclines,0.,fit)
    for i in range(w.size):
        diff = abs(w[i]-modellines)
        if diff.min()<5.*scale:
            left_matches.append([arclines[i],modellines[diff.argmin()]])

    fit = sf.lsqfit(scipy.asarray(left_matches),'polynomial',order)


    """ Optimization function for refining the wavelength solution. """
    def dofit(p,x,data,model):
        fit = {'coeff':scipy.atleast_2d(p).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        m = interpolate.splev(w,model)
        return (m-data)

    x = scipy.arange(arc.size).astype(scipy.float32)

    initial_pars = fit['coeff'][:,0].tolist()
    coeff,ier = optimize.leastsq(dofit,initial_pars,
                        (x,arc_wide,arcmodel['wide']),maxfev=100000)
    coeff,ier = optimize.leastsq(dofit,coeff,
                        (x,arc,arcmodel['matched']),maxfev=100000)
    outcoeff = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}

    def skycorrect(p,arc,sky,arcmodel,skymodel):
        fit = {'coeff':scipy.atleast_2d(p[:-1]).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        arcm = interpolate.splev(w+p[-1],arcmodel)
        chi_arc = (arcm-arc)
        s = sky[w>5100.]
        skym = interpolate.splev(w[w>5100.],skymodel)
        skym *= scipy.median(s/skym)
        chi_sky = 5.*(skym-s)#/abs(m)**0.5
        chi = scipy.concatenate((chi_arc,chi_sky))
        return chi

    newcoeff = coeff.tolist()
    newcoeff.append(0.)
    coeff,ier = optimize.leastsq(skycorrect,newcoeff,
                        (arc,sky,arcmodel['matched'],skymodel['matched']),
                        maxfev=100000)
    outcoeff = {'coeff':scipy.atleast_2d(coeff[:-1]).T,'type':'polynomial'}

    """
    wave = sf.genfunc(x,0.,outcoeff)
    sky = sky[wave>5000.]
    wave = wave[wave>5000.]

    m = interpolate.splev(wave,wavemodel['matched'])
    ratio = scipy.median(sky/m)
    import pylab
    pylab.plot(wave,sky)
    pylab.plot(wave,m*ratio)
    pylab.show()

    offset,ier = optimize.leastsq(skycorrect,[0.],
                        (wave,sky,wavemodel['matched']),maxfev=100000)
    print offset
    outcoeff['coeff'][0] += offset
    """

    return outcoeff


def wave_skylines(sky,solution):
    STD_LINES = [5197.928,5200.286,5202.977,
                 5460.735,5577.345,5867.5522,5915.308,5932.864,6257.970,
                 6300.320,6363.810,
                 6533.040,6553.610,6863.971,6912.620,6923.210,6939.520,
                 7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,
                 7821.510,7841.266,7993.332,8310.719,8344.613,8399.160,8415.231,
                 8430.170,8791.186,8885.830,8943.395,8988.384,9038.059,9337.854,
                 9375.977,9419.746,9439.670,9458.524]
    x = scipy.arange(sky.size).astype(scipy.float32)
    lines = get_lines(x,sky)

    scale = solution['coeff'][1]
    order = solution['coeff'].size - 1

    if scale>1.5:
       STD_LINES.insert(6,5891.)

    w = sf.genfunc(lines,0.,solution)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    lines = get_lines(x,sky,nstd=7.)
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    return fit


def wave_arclines(arc,arcmodel,sky,solution):
    from scipy import interpolate
    STD_LINES = scipy.sort(arcmodel['lines'])
    SKYLINES = [5577.338,6300.304,6363.78,6553.617]#,6912.62]
    x = scipy.arange(arc.size).astype(scipy.float32)
    lines = get_lines(x,arc)
    skycoords = get_lines(x,sky)

    arccoords = lines.copy()
    scale = arcmodel['scale']

    if scale>1.5:
        SKYLINES.insert(1,5891.)

    order = solution['coeff'].size - 1
    fit = solution
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for line in STD_LINES:
        diff = abs(w-line)
        if diff.min()<5.*scale:
            pos = diff.argmin()
            matches.append([lines[pos],line])
    print matches,order
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    order = 1

    for i in range(7):
        matched = [m for m in matches]
        if skycoords.size==0:
            offset = 0.
            break
        skyline = sf.genfunc(skycoords,0.,fit)
        offsets = []
        for line in SKYLINES[:i*2+2]:
            diff = abs(line-skyline)
            if diff.min()<5.*scale:
                offsets.append(line-skyline[diff.argmin()])
        if len(offsets)==0:
            offset = 0.
            break
        offset = scipy.median(scipy.asarray(offsets))

        for line in SKYLINES[:i*2+2]:
            diff = abs(line-skyline)
            if diff.min()<5.*scale:
                pos = diff.argmin()
                matched.append([skycoords[pos],line-offset])
        fit = sf.lsqfit(scipy.asarray(matched),'polynomial',order)


    fit['coeff'][0] += offset

    return fit


def wave_arcsky(arc,arcmodel,sky,solution):
    """
    First find the best solution with the skylines, then apply this solution
        to all arclines within the bounds of the lowest/highest wavelength
        skylines, solving for the delta_pixel offset between the sky and the
        arc. Then find the solution for all lines (sky and delta_pixel-offset
        arcs).
    """
    def clip(arr):
        a = arr.copy()
        m,s,l = a.mean(),a.std(),a.size
        while 1:
            a = a[abs(a-m)<3.*s]
            if a.size==l:
                return m,s
            m,s,l = a.mean(),a.std(),a.size

    STD_LINES = [5197.928,5200.286,5202.977,
                 5460.735,5577.345,5867.5522,5915.308,5932.864,6257.970,
                 6300.320,6363.810,
                 6533.040,6553.610,6863.971,6912.620,6923.210,6939.520,
                 7303.716,7329.148,7340.885,7358.659,7392.198,7586.093,7808.467,
                 7821.510,7841.266,7993.332,8310.719,8344.613,8399.160,8415.231,
                 8430.170,8791.186,8885.830,8943.395,8988.384,9038.059,9337.854,
                 9375.977,9419.746,9439.670,9458.524]
    x = scipy.arange(sky.size).astype(scipy.float32)
    lines = get_lines(x,sky)

    global fit
    scale = solution['coeff'][1]
    order = solution['coeff'].size - 1

    if scale>1.5:
        STD_LINES.insert(6,5891.)

    w = sf.genfunc(lines,0.,solution)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<5.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(matches),'polynomial',order)

    lines = get_lines(x,sky,nstd=7.)
    w = sf.genfunc(lines,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    matches = scipy.asarray(matches)
    fit = sf.lsqfit(matches,'polynomial',order)
    revfit = sf.lsqfit(matches[:,::-1],'polynomial',order)

    ARCS = scipy.sort(arcmodel['lines'])
    alines = get_lines(x,arc)
    xmin,xmax = matches[0,0],matches[-1,0]
    arc_x = sf.genfunc(ARCS,0.,revfit)
    offset = []
    for i in range(arc_x.size):
        if arc_x[i]<xmin-2 or arc_x[i]>xmax+2:
            continue
        diff = arc_x[i]-alines
        if abs(diff).min()<9.:
            offset.append(diff[abs(diff).argmin()])
    offset = scipy.asarray(offset)
    off,width = clip(offset)

    aw = sf.genfunc(alines+off,0.,fit)
    matches = []
    for i in range(w.size):
        diff = abs(w[i]-STD_LINES)
        if diff.min()<3.*scale:
            matches.append([lines[i],STD_LINES[diff.argmin()]])
    k = len(matches)
    for i in range(aw.size):
        diff = abs(aw[i]-ARCS)
        if diff.min()<3.*scale:
            matches.append([alines[i]+off,ARCS[diff.argmin()]])
    matches = scipy.asarray(matches)
    fit = sf.lsqfit(matches,'polynomial',order)
    """
    # This code is to optimize the arc offset -- likely unnecessary.
    def opt(p):
        global fit
        off = p[0]
        w = sf.genfunc(lines,0.,fit)
        aw = sf.genfunc(alines+off,0.,fit)
        matches = []
        for i in range(w.size):
            diff = abs(w[i]-STD_LINES)
            if diff.min()<3.*scale:
                matches.append([lines[i],STD_LINES[diff.argmin()]])
        k = len(matches)
        for i in range(aw.size):
            diff = abs(aw[i]-ARCS)
            if diff.min()<3.*scale:
                matches.append([alines[i]+off,ARCS[diff.argmin()]])
        matches = scipy.asarray(matches)
        fit = sf.lsqfit(matches,'polynomial',order)
        return (matches[:,1]-sf.genfunc(matches[:,0],0.,fit))
    from scipy import optimize
    coeff,ier = optimize.leastsq(opt,[off],epsfcn=1e-15)
    """
    return fit


def get_lines(xvals,data,dosat=True,std=None,nstd=15.):
    from scipy import ndimage
    def clip(arr,nsig):
        a = arr.copy()
        m,std,size = a.mean(),a.std(),a.size

        while 1:
            flag = a*0.
            cond = abs(a-m)>nsig*std
            flag[cond] = 1.
            flag = ndimage.maximum_filter(flag,5)
            a = a[flag==0]
            if a.size==size or a.size<2:
                return m,std
            m,std,size = a.mean(),a.std(),a.size

    bg = ndimage.percentile_filter(data,10,51)

    if std is None:
        m,std = clip(data-bg,3.)
    maxfilt = ndimage.maximum_filter(data,9)
    # nstd: how many std above bg should we identify as the trace
    # try a lower value of nstd if we don't have a bright star trace...
    peaks = scipy.where((maxfilt==data)&((data-bg)>nstd*std))[0]

    out = []
    fit = scipy.empty(4)
    for p in peaks:
        if p-14<0 or p+15>xvals.size:
            continue
        fit[0] = bg[p]
        fit[1] = data[p]-fit[0]
        fit[2] = xvals[p]
        fit[3] = 1.

        if data[p]<54000.:
            fitdata = scipy.empty((9,2))
            fitdata[:,0] = xvals[p-4:p+5].copy()
            fitdata[:,1] = data[p-4:p+5].copy()-bg[p-4:p+5].mean()
        elif dosat:
            fitdata = scipy.empty((25,2))
            fitdata[:,0] = xvals[p-12:p+13].copy()
            fitdata[:,1] = data[p-12:p+13].copy()-bg[p-12:p+13]
            fitdata = fitdata[fitdata[:,1]<54000.]
        else:
            continue

        answer = (fitdata[:,0]*fitdata[:,1]).sum()/fitdata[:,1].sum()
        out.append(answer)

    return scipy.asarray(out)


def combine_xw(coords,xsoln,wsoln,xord,yord):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newx = sf.genfunc(x,y,xsoln['back'])
    wave = sf.genfunc(newx,0.,wsoln)

    data = scipy.empty((x.size,3))
    data[:,0] = x.copy()
    data[:,1] = y.copy()
    data[:,2] = wave.copy()
    soln = sf.lsqfit(data,'chebyshev',xord,yord)

    tmp = data[:,0].copy()
    data[:,0] = data[:,2].copy()
    data[:,2] = tmp.copy()
    soln2 = sf.lsqfit(data,'chebyshev',xord,yord)

    return {'back':soln,'forw':soln2}


def combine_xy(offset,coords,xsoln,ysoln,xord,yord):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newy = ysoln['back'].flatten()
    newx = sf.genfunc(x,y-offset,xsoln['back'])

    data = scipy.empty((x.size,3))
    data[:,0] = x.copy()
    data[:,1] = y.copy()
    data[:,2] = newx.copy()
    soln = sf.lsqfit(data,'chebyshev',xord,yord)

    tmp = data[:,0].copy()
    data[:,0] = data[:,2].copy()
    data[:,2] = tmp.copy()
    soln2 = sf.lsqfit(data,'chebyshev',xord,yord)

    return {'back':soln,'forw':soln2}


def combine_xyw(coords,xsoln,ysoln,wsoln,xord,yord):
    x = coords[1].flatten()
    y = coords[0].flatten()

    newy = ysoln.flatten()
    from scipy import random
    k = random.random(x.size)
    args = k.argsort()
    x = x[args[:y.size/10]]
    newy = newy[args[:y.size/10]]
    y = y[args[:y.size/10]]

    newx = sf.genfunc(x,y,xsoln['back'])
    wave = sf.genfunc(newx,0.,wsoln)

    data = scipy.empty((x.size,3))
    data[:,0] = wave.copy()
    data[:,1] = y.copy()
    data[:,2] = x.copy()
    output_to_ccdx = sf.lsqfit(data,'chebyshev',xord,yord)

    data[:,2] = newy.copy()
    output_to_ccdy = sf.lsqfit(data,'chebyshev',xord,yord)

    data[:,0] = x.copy()
    data[:,1] = y.copy()
    data[:,2] = wave.copy()
    ccdx_ycor_to_wave = sf.lsqfit(data,'chebyshev',xord,yord)

    return {'sky2x':output_to_ccdx,'sky2y':output_to_ccdy,'ccd2wave':ccdx_ycor_to_wave}


def arcwave2(arc,arcmodel,scale,order=3,bcutoff=2e3,rcutoff=1e4):
    from scipy import ndimage,stats,interpolate,signal,optimize
    import pylab

    arc = arc.copy()
    arc[scipy.isnan(arc)] = 1.
    arc[arc<=1.] = arc.mean()

    wave,model = arcmodel['orig']
    model[scipy.isnan(model)] = 0.
    wave = wave.copy()
    model = model.copy()
    cond = (wave>bcutoff)&(wave<rcutoff)
    corrmodel = model.copy()
    corrmodel[(wave<bcutoff)|(wave>rcutoff)] = 0.

    corr = signal.correlate(arc,corrmodel,mode='valid')
    offset = corr.argmax()

    lines = arcmodel['lines'].copy()
    bc = lines.min()-scale*20.
    rc = lines.max()+scale*20.

    x = scipy.arange(arc[offset:].size)
    w = wave[offset:offset+x.size].copy()
    cond = (w>bc)&(w<rc)
    fit = scipy.empty((x[cond].size,2))
    fit[:,0] = x[cond].copy()
    fit[:,1] = w[cond].copy()
    pars = sf.lsqfit(fit,'polynomial',1)

    pars['coeff'][0] = wave[offset]
    pars['coeff'][1] = scale
#    pars = [wave[offset],scale]
#    for i in range(2,order+1):
#        pars.append(1e-5**i)
    pylab.plot(wave,model)
    w = sf.genfunc(scipy.arange(x.size),0.,pars)
    pylab.plot(w,interpolate.splev(w,arcmodel['matched']))
    pylab.show()

    def arcfit(p,x,arc,mod):
        fit = {'coeff':scipy.atleast_2d(p).T,'type':'polynomial'}
        w = sf.genfunc(x,0.,fit)
        cond = (w>bcutoff)&(w<rcutoff)
        m = interpolate.splev(w[cond],mod)
        chi = (m-arc[cond])/abs(arc[cond])**0.5
        return chi

    widearc = ndimage.gaussian_filter(arc,7.)
    x = scipy.arange(arc[offset:].size)
    coeff,ier = optimize.leastsq(arcfit,pars,(x,widearc[offset:].copy(),
                        arcmodel['wide']),maxfev=100000)

    fit = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}


    x = scipy.arange(arc.size)
    l = get_lines(x,arc,nstd=15.)
    lw = sf.genfunc(l-offset,0.,fit)
    lines = []
    for i in range(l.size):
        diff = abs(lw[i]-arcmodel['lines'])
        if diff.min()>5.*scale:
            continue
        lines.append([l[i],arcmodel['lines'][diff.argmin()]])
    fit = sf.lsqfit(scipy.asarray(lines),'polynomial',order)

    pars = fit['coeff'].flatten()
    coeff,ier = optimize.leastsq(arcfit,pars,(x[offset:],arc[offset:],
                        arcmodel['matched']),maxfev=100000)

    fit = {'coeff':scipy.atleast_2d(coeff).T,'type':'polynomial'}

    return fit
