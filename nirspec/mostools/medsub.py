import pyfits,glob,numpy
import special_functions as sf
from scipy import ndimage,stats

def medsubtract(image,outname):
    data = pyfits.open(image)[0].data.copy()
    if data.ndim==3:
        data = data[0].copy()

    tmp = data.copy()
    tmp[numpy.isnan(tmp)] = 0.
    tmp -= numpy.sort(tmp,0)[tmp.shape[0]/5]
    trace = tmp.sum(1)
    peak = trace.argmax()

    center = numpy.empty(data.shape[1])
    w = center.copy()
    for i in range(1+center.size/100):
       b = i*100
       e = b+100
       if e>center.size:
           e = center.size
       if b==e:
           continue
       center[b:e] = tmp[:,b:e].sum(1).argmax()
    bg = center.copy()
    x = numpy.arange(data.shape[0])
    for i in range(center.size):
        d = tmp[:,i].copy()
        peak = center[i]
        if numpy.isnan(d[peak]):
            center[i] = peak
            continue
        fit = numpy.array([0.,d[peak],peak,1.])
        cond = ~numpy.isnan(d)
        input = numpy.empty((d[cond].size,2))
        input[:,0] = x[cond].copy()
        input[:,1] = d[cond].copy()
        fit,chi = sf.ngaussfit(input,fit)
        center[i] = fit[2]
        w[i] = fit[3]

    fit = sf.lsqfit(ndimage.median_filter(center,17),'polynomial',5)
    centroid = sf.genfunc(numpy.arange(bg.size),0.,fit)
    w = numpy.median(w)
    for i in range(bg.size):
        d = data[:,i].copy()
        d[centroid[i]-w*4:centroid[i]+w*4] = numpy.nan
        data[:,i] -= stats.nanmedian(d)

    hdu = pyfits.open(image)[0]
    hdu.data = data.copy()
    hdu.writeto(outname,clobber=True)
