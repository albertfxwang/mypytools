import pyfits,pylab,scipy,os,sys,glob,numpy,cPickle
from mostools import spectools as st
import special_functions as sf
from scipy import ndimage,interpolate


def extract(image,varimage,outname,width=2.,offset=0.,pos=None,minwave=None,maxwave=None,response_image=None,response_model=None,regions=[],centroid=None):

    if outname.lower().find('fits')<1:
        outname = "%s.fits"%outname

    resp = None
    respwave = None
    if response_image is not None:
        if response_model is None:
#            print "No response model included; not performing response correction."
            import cPickle
            resp = cPickle.load(open(response_image))
        else:
            resp = get_response(response_image,response_model,regions)
    elif response_model is not None:
        print "No response image included; not performing response correction."


    if minwave is None:
        minwave = -10.
    if maxwave is None:
        maxwave = 1e99

    data = pyfits.open(image)[0].data.copy()
    if data.ndim==3:
        if data.shape[0]>1:
            print "Only working on first image plane of 3D image."
        data = data[0].copy()
    wave = st.wavelength(image)
    indx = scipy.where((wave>minwave)&(wave<maxwave))[0]
    minx,maxx = indx.min(),indx.max()
    wave = wave[minx:maxx]
    data = data[:,minx:maxx]
    tmp = data.copy()

    tmp[numpy.isnan(data)] = 0.

    ospec = wave*0.
    ovar = wave*0.
    x = numpy.arange(data.shape[0])
    if centroid is None:
        n = tmp.shape[1]/100
        peaks = []
        bounds = numpy.linspace(0,tmp.shape[1],n+1)
        for i in range(n):
            a,b = bounds[i],bounds[i+1]
            p = tmp[:,a:b].sum(1).argmax()
            peaks.append([(a+b)/2.,p])
        peak = sf.lsqfit(numpy.asarray(peaks),'polynomial',3)

        if pos is not None:
            peak['coeff'][0] = pos

        peaks = sf.genfunc(numpy.arange(tmp.shape[1]),0.,peak)

        center = wave*0.
        for i in range(center.size):
            d = data[:,i].copy()
            peak = peaks[i]
            if peak<0:
                peak = 0
            if peak>=d.shape:
                peak = d.shape-1.
            fit = numpy.array([0.,d[peak],peak,1.])
            cond = ~numpy.isnan(d)
            input = numpy.empty((d[cond].size,2))
            input[:,0] = x[cond].copy()
            input[:,1] = d[cond].copy()
            fit = sf.ngaussfit(input,fit)[0]
            center[i] = fit[2]

        fit = sf.lsqfit(ndimage.median_filter(center,17),'polynomial',5)
        centroid = sf.genfunc(scipy.arange(wave.size),0.,fit)
    #import pylab
    #pylab.plot(center-centroid)
    #pylab.show()

    xvals = scipy.arange(data.shape[0])
    var = pyfits.open(varimage)[0].data.copy()
    if var.ndim==3:
        var = var[0].copy()
    var = var[:,minx:maxx]
    cond = (numpy.isnan(var))|(numpy.isnan(data))

    for i in range(ovar.size):
        c = centroid[i]+offset
        wid = width

        mask = xvals*0.
        mask[(xvals-c<wid)&(xvals-c>0.)] = wid-(xvals[(xvals-c<wid)&(xvals-c>0.)]-c)
        mask[(xvals-c<wid-1)&(xvals-c>0.)] = 1.
        mask[(c-xvals<wid+1)&(c-xvals>0.)] = (wid+1)-(c-xvals[(c-xvals<wid+1)&(c-xvals>0.)])
        mask[(c-xvals<wid)&(c-xvals>0.)] = 1.
        mask[cond[:,i]] = 0.
        mask /= mask.sum()

        ospec[i] = (mask[~cond[:,i]]*data[:,i][~cond[:,i]]).sum()
        ovar[i] = ((mask[~cond[:,i]]**2)*var[:,i][~cond[:,i]]).sum()

    badpix = numpy.where(numpy.isnan(ospec))[0]
    ospec[badpix] = 0.
    ovar[badpix] = ovar[~numpy.isnan(ovar)].max()*1e6

    med = numpy.median(ospec)
    ospec /= med
    ovar /= med**2

    if resp is not None:
#        if respwave is None:
#            resp = sf.genfunc(wave,0,resp)
#        else:
#        resp = interpolate.splrep(respwave,resp)
        resp = interpolate.splev(wave,resp)
        ospec *= resp
        ovar *= resp**2

    st.make_spec(ospec,ovar,wave,outname,clobber=True)
    return centroid

def get_response(image,model,regions):
    model = numpy.loadtxt(model)
    mwave = model[:,0]
    mspec = model[:,1]

    model = interpolate.splrep(mwave,mspec,k=3,s=0)

    hdu = pyfits.open(image)
    if len(hdu)==4:
        spec = hdu[1].data.copy()
        wave = st.wavelength(image,1)
    else:
        spec = hdu[0].data.copy()
        wave = st.wavelength(image)

    outmodel = interpolate.splev(wave,model)

    ratio = outmodel/spec

    badregions = []
    cond = ~numpy.isnan(ratio)
    cond = cond&(~numpy.isinf(ratio))
    for lo,hi in regions:
        badregions.append([lo,hi])
        cond = cond&(~((wave>lo)&(wave<hi)))
    scurrent = 2.*wave[cond].size**0.5
    smod = ratio[cond].mean()**2
    spmodel = interpolate.splrep(wave[cond],ratio[cond],k=3,s=scurrent*smod)
    resp = None
    while resp!='q' and resp!='Q':
        import pylab
        current = interpolate.splev(wave,spmodel)
        pylab.plot(wave,current)
        pylab.plot(wave,ratio)
        pylab.gca().fmt_xdata = pylab.FormatStrFormatter('%7.2f')
        pylab.show()
        resp = raw_input("Enter command (q, m, s, w, h): ")
        if resp=='m':
            region = raw_input("Enter region to mask (eg, 6530,6580): ").split(',')
            while len(region)!=2:
                region = raw_input("Please input wavelengths joined by a comma: ").split(',')
            lo,hi = float(region[0]),float(region[1])
            badregions.append([lo,hi])
            cond = cond&(~((wave>lo)&(wave<hi)))
            spmodel = interpolate.splrep(wave[cond],ratio[cond],k=3,s=scurrent*smod)
        elif resp=='s':
            scurrent = float(raw_input("Current smoothing factor is %4.2f, enter new smoothing factor: "%scurrent))
            spmodel = interpolate.splrep(wave[cond],ratio[cond],k=3,s=scurrent*smod)
        elif resp=='h':
            print "Use q to quit, m to mask, s to set smoothing scale, w to write model to disk"
        elif resp=='w':
            import cPickle
            outname = raw_input("Name of file to write to: ")
            f = open(outname,'wb')
            cPickle.dump(spmodel,f)
            f.close()

    print "Regions masked:",badregions
    return spmodel
