import pylab,scipy,pyfits,special_functions,numpy
import matplotlib as mpl

class SpecPlot:

    def __init__(self,file):
        f = pyfits.open(file)
        self.filename = file

        """
        Test to see if data are in MWA format (ie ext0=meta, ext1=sci,
          ext2=smooth, ext3=var). If not, set smooth and variance data
          to science data.
        """
        self.ismwa = True
        self.isSDSS = False
        self.isSpec1d = False
        if len(f)==5:
            self.ismwa = False
            self.isSpec1d = True
        else:
            try:
                self.position = f[0].header['center']
            except:
                self.ismwa = False
                try:
                    nspex = f[0].data.shape[0]
                    if nspex==5:
                        self.isSDSS = True
                except:
                    pass
        if self.ismwa:
            self.pos = self.position
            self.width = f[0].header['width']
            self.rawdata = f[1].data.astype(scipy.float64)
            self.current = f[2].data.astype(scipy.float64)
            self.varspec = f[3].data.astype(scipy.float64)
            f.close()
            self.wave = wavelength(file,1)
        elif self.isSDSS:
            self.rawdata = f[0].data[0].copy()
            self.current = f[0].data[0].copy()
            self.varspec = f[0].data[2]**2
            self.wave = wavelength(file)
        elif self.isSpec1d:
            self.rawdata = numpy.concatenate((f[1].data.SPEC[0],f[2].data.SPEC[0]))
            self.current = numpy.concatenate((f[1].data.SPEC[0],f[2].data.SPEC[0]))
            self.varspec = numpy.concatenate((f[1].data.IVAR[0]**-1,f[2].data.IVAR[0]**-1))
            self.wave = numpy.concatenate((f[1].data.LAMBDA[0],f[2].data.LAMBDA[0]))

        else:
            self.rawdata = f[0].data.astype(scipy.float64)
            self.current = f[0].data.astype(scipy.float64)
            self.varspec = f[0].data.astype(scipy.float64)
            f.close()
            self.wave = wavelength(file)

            self.rawdata = self.rawdata.reshape(self.rawdata.size)
            self.current = self.current.reshape(self.current.size)
            self.varspec = self.varspec.reshape(self.varspec.size)

        """
        Trim bad edges.
        """
        tmp = scipy.where((self.rawdata!=0)&~scipy.isnan(self.rawdata))[0]
        left = tmp.min()
        right = tmp.max()+1
        self.current = self.current[left:right]
        self.varspec = self.varspec[left:right]
        self.wave = self.wave[left:right]
        self.rawdata = self.rawdata[left:right]


        """
        Linelist for galaxies.
        """
        self.linelist = [2344,2374,2383,2587,2599.4,2750.3,2795.,2802.7,2852,3726.03,3728.82,3835.38,3889.05,3933.66,3968.47,4101.74,4305.,4340.47,4861.33,4962.,5007.,5177.,5270.,5328.,5341.,5371.5,5891.9,6562.8,6583.5,6716.4,6730.8,8498.,8542.,8662.]
        self.linename = ['Fe','','Fe','Fe','','FeII','','MgII','MgI',"[OII]","",'Heta','','CaK','CaH','H-delta','Gband','H-gamma','H-beta','OIII','OIII','Mgb','Fe','Fe','','Fe','NaII','H-alpha','N','S','','CaT','CaT','CaT']

        """
        Linelist for high-z quasars.
        """
        self.qso = False
        self.qsolist = [1215.24,1239.4,1305.53,1335.52,1397.61,1399.8,1545.86,1637.85,1665.85,1857.4,1908.27,2326.0,2439.5,2800.32]
        self.qsoname = ['Lya','Nv','Oi','Cii','SiIV','Oiv','Civ','HeII','Oiii','AlIII','Ciii','Cii','NeIV','MgII']

        self.collect = [] # Stores ID'd lines

        self.skyid = None
        self.snid = None

        self.z = 0.

        self.axes = None
        self.redraw()
        self.axes.fmt_xdata = pylab.FormatStrFormatter("%4.2f")

    def redshift(self, z):
        # Setup transformations stuff
        #blend = mpl.transforms.blend_xy_sep_transform
        blend = mpl.transforms.blended_transform_factory
        trans = blend(self.axes.transData,self.axes.transAxes)
        lines = self.linelist
        names = self.linename

        if self.qso:
            lines = self.qsolist
            names = self.qsoname

        self.axes.texts = []
        self.axes.collections = []
        axis = self.axes.axis()
        ymin,ymax = axis[2],axis[3]
        if self.skyid==None:
            pass
        else:
            sky = self.axes.lines[self.skyid]
        self.axes.lines = [self.axes.lines[0]]
        if self.skyid==None:
            pass
        else:
            self.axes.lines.append(sky)
            self.skyid = 1
        plotwaves = []
        for i in range(len(lines)):
            wave = lines[i]
            wave *= (1.+z)
            if wave>self.wave[-1] or wave<self.wave[0]:
                continue
            plotwaves.append(wave)
            pylab.axvline(wave,color='black')
            self.axes.text(wave,1,names[i],transform=trans,rotation=60,ha='center')
#            pylab.text(wave,axis[3],names[i],rotation=60)
#        if len(plotwaves)>0:
#            self.axes.vlines(plotwaves,ymin,ymax,color='black')
        lines = [3797,3868,4383,4455,4531]
        plotwaves = []
        for i in range(len(lines)):
            wave = lines[i]*(1.+z)
            if wave>self.wave[-1] or wave<self.wave[0]:
                continue
            plotwaves.append(wave)
            pylab.axvline(wave,color='red')
#        if len(plotwaves)>0:
#            self.axes.vlines(plotwaves,ymin,ymax,color='red')
#            self.axes.axvline(plotwaves,0,1,color='red')
        self.z = z
        pylab.draw()

    def id(self,wave,feature):
        # Get closest feature to input value
        lines = scipy.asarray(self.linelist)
        diff = abs(lines-feature)
        if diff.min()<3.:
            feature = lines[diff.argmin()]

        disp = self.wave[1]-self.wave[0]
        indx = scipy.where(abs(wave-self.wave)<disp)[0][0]
        if indx<5 or indx>self.wave.size-4:
            self.lineid(wave,feature)
            return
        fitdata = scipy.empty((11,2))
        fitdata[:,0] = self.wave[indx-5:indx+6]
        fitdata[:,1] = self.current[indx-5:indx+6]

        fit = scipy.zeros(4)
        fit[1] = self.current[indx]
        fit[2] = wave
        fit[3] = 2.

        fit,chi2 = special_functions.ngaussfit(fitdata,fit)
        if abs(fit[2]-wave)>5.:
            self.lineid(wave,feature)
        else:
            self.lineid(fit[2],feature)
        self.collect.append(self.z)

    def erase_z(self):
        self.collect = []

    def get_z(self):
        if len(self.collect)==0:
            return self.z
        return scipy.asarray(self.collect).mean()

    def get_z_err(self):
        if len(self.collect)<2:
            return 0.
        return scipy.asarray(self.collect).std()

    def lineid(self,wave,feature):
        z = (float(wave)/feature)-1.
        print "wave = %7.2f, z = %2.6f" % (wave,z)
        self.redshift(z)

    def smooth(self,width=5):
        from scipy import signal
        var = self.varspec.copy()
        var[var<=0] = 1.
        var[scipy.isnan(var)] = 1.
        self.current = signal.wiener(self.current,width,var)
        self.axes.lines[0].set_ydata(self.current)
        self.axes.set_ylim([self.current.min(),self.current.max()])
        pylab.draw()

    def gauss(self,width=5):
        from scipy import ndimage
        self.current = ndimage.gaussian_filter(self.current,width)
        self.axes.lines[0].set_ydata(self.current)
        self.axes.set_ylim([self.current.min(),self.current.max()])
        pylab.draw()

    def boxcar(self,width=5):
        from scipy import ndimage
        self.current = ndimage.uniform_filter(self.current,width)
        self.axes.lines[0].set_ydata(self.current)
        self.axes.set_ylim([self.current.min(),self.current.max()])
        pylab.draw()

    def close(self):
        pylab.close()

    def plot(self):
        pylab.close()
        pylab.plot(self.wave,self.current,linestyle='steps')
        pylab.show()
        self.axes = pylab.gca()
        self.axes.fmt_xdata = pylab.FormatStrFormatter("%4.2f")

    def sky(self):
        if self.skyid is None:
            skymax = scipy.nanmax(self.varspec)
            scimax = scipy.nanmax(self.current)
            scimin = scipy.nanmin(self.current)
            spec = self.varspec*scimax/skymax - abs(scimin)
            self.skyid = len(self.axes.lines)
            pylab.plot(self.wave,spec,color='green')
        else:
            self.axes.lines.remove(self.axes.lines[self.skyid])
            self.skyid = None
            pylab.draw()

    def sn(self):
        if self.snid is None:
            self.snid = len(self.axes.lines)
            pylab.plot(self.wave,self.rawdata/self.varspec**0.5)
        else:
            self.axes.lines.remove(self.axes.lines[self.snid])
            self.snid = None
            pylab.draw()

    def redraw(self):
        self.current = self.rawdata.copy()
        self.plot()
        self.axes = pylab.gca()
        self.skyid = None

    def clear(self):
        self.axes.texts = []
        self.axes.lines = [self.axes.lines[0]]
        pylab.draw()
        self.skyid = None

def parse_hdr(header):
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
