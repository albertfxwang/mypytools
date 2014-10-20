import numpy,pyfits,cPickle,sys
import special_functions as sf
from mostools import spectools as st
from scipy import optimize,ndimage,interpolate

dir = st.__file__.split('spectools.py')[0]

def vegacorr(input,output,size=250):
    d = pyfits.open(input)[1].data.copy()
    wave = st.wavelength(input,1)
    k = numpy.linspace(wave[0],wave[-1],wave.size/size)

    twave,trans = cPickle.load(open(dir+'data/trans.dat'))
    kernel = (wave[1]-wave[0])/(1e4*(twave[1]-twave[0]))
    trans = ndimage.gaussian_filter(trans,kernel)
    sky = interpolate.splrep(twave*1e4,trans,k=3,s=0)

    swave,star = cPickle.load(open(dir+'data/vega.dat'))
    snorm = star.mean()
    starmodel = interpolate.splrep(swave*10.,star/snorm,k=3,s=0)

    def dofit(p,wave,star,vega,sky,trans=False):
        vvel,broad,svel,scale = p
        vwave = 10**(numpy.log10(wave)*vvel)
        v = interpolate.splev(vwave,vega)
        v = ndimage.gaussian_filter(v,broad)
        swave = 10**(numpy.log10(wave)*svel)
        s = interpolate.splev(swave,sky)**scale
        s[swave<8510.] = 1.
        mod = interpolate.splrep(wave,star/(v*s),t=k[1:-1],k=3,task=-1)
        model = s*v*interpolate.splev(wave,mod)
        if trans:
            return v/model
        return model-star

    pars = [1.,1.,1.,1.]
    coeff,ier = optimize.leastsq(dofit,pars,(wave,d,starmodel,sky))

    model = dofit(coeff,wave,d,starmodel,sky,True)/snorm
    model = interpolate.splrep(wave,model)
    f = open(output,'wb')
    cPickle.dump(model,f,2)
    f.close()
    return model
