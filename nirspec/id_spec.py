import scipy,pickle,numpy
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    import pylab as plt
import special_functions as sf
from scipy import ndimage
STANDARD = None
print "plt.__file__:", plt.__file__

class IDSpectrum:
    """
    IDSpectrum class for identification of spectral lines starting with an
        initial model of the spectrum.
    """
    def __init__(self, data,sky,wave,wave2,skymodel):
        """ Plot data """
        self.data = plt.plot(wave,data,c='b')[0]
        self.sky = plt.plot(wave2,sky,c='gray')[0]
        self.canvas = self.data.get_figure().canvas
        self.ax = self.data.get_axes()
        self.xdata = self.data.get_xdata().copy()
        self.start = [data.copy(),sky.copy(),wave.copy(),wave2.copy(),skymodel]
        self.skymodel = skymodel

        """ Get metadata (ie line locations) for arcs """
        data = self.data.get_ydata()
        self.datapeaks = ndimage.maximum_filter(data,9)
        tmp = scipy.sort(data)
        thresh = tmp[tmp.size*0.95]
        cond = (data==self.datapeaks)&(data>thresh)
        self.datapeaks = scipy.where(cond)[0]
        self.datasel = self.datapeaks*0
        self.datalines = []
        for peak in self.datapeaks:
            l = plt.axvline(self.xdata[peak],c='k')
            self.datalines.append(l)

        self.spec = self.data
        self.peaks = self.datapeaks
        self.selected = self.datasel
        self.lines = self.datalines

        """ Set useful flags """
        self.domotion = False
        self.origx = None
        self.soln = None
        self.pick = False
        self.fitlines = None
        self.keyid = self.canvas.mpl_connect('key_press_event',self.key_press)

        self.connect()

        print """
Mouse Controls:
    - left button drags single lines (rescales spectrum!)
    - middle button drags all lines (or exits from pan/zoom modes)
    - right button selects/deselects lines
Keyboard Commands:
    a - add new line (use mouse to select the line)
    m - fit a polynomial to the blue `solution'
    d - optimize the blue fit to the gray model (like m, but optimizes too)
    w - write the current state to disk
    r - read a saved state
    n - reset to the initial state
    q - quit (performs an `m' fit if no fit has been applied yet)
"""
        plt.show()

    def connect(self):
        """ Connect the mouse to the plot """
        self.pressid = self.canvas.mpl_connect('button_press_event',
                                                self.on_press)
        self.moveid = self.canvas.mpl_connect('motion_notify_event',
                                                self.on_motion)
        self.offid = self.canvas.mpl_connect('button_release_event',
                                                self.on_release)

    def on_press(self,event):
        """
        Deal with mouse button presses, including stretching, shifting,
            and line identification.
        """
        """ Turn off plot tools """
        if self.canvas.toolbar.mode!='':
            if event.button==2:
                self.canvas.toolbar.zoom()
                self.canvas.toolbar.pan()
                self.canvas.toolbar.pan()
            return

        self.xdata = self.spec.get_xdata().copy()
        if event.xdata==None or event.ydata==None: # Not in axes
            return

        ind = abs(self.xdata-event.xdata).argmin()
        indx = abs(self.peaks-ind).argmin()

        if abs(self.peaks-ind).min()>4.: # Not near any peaks
            return

        """ Select/unselect lines """
        if event.button==3:
            if self.selected[indx]==1:
                self.selected[indx] = 0
                self.lines[indx].set_color('k')
            else:
                self.selected[indx] = 1
                self.lines[indx].set_color('r')
            plt.draw()
            return


        self.origx = self.xdata[self.peaks[indx]]
        self.lbound = self.xdata[0]
        self.lindx = 0
        for i in range(indx):
            if self.xdata[self.peaks[i]]>self.origx:
                break
            if self.selected[i]==1:
                self.lbound = self.xdata[self.peaks[i]]
                self.lindx = self.peaks[i]

        self.rbound = self.xdata[-1]
        self.rindx = -1
        for i in range(indx+1,self.peaks.size):
            if self.xdata[self.peaks[i]]<self.origx:
                continue
            if self.selected[i]==1:
                self.rbound = self.xdata[self.peaks[i]]
                self.rindx = self.peaks[i]
                break
        self.rbound = 1e6
        self.lbound = 1
        self.selected[indx] = 1
        self.indx = self.peaks[indx]
        self.domotion = True
        self.lines[indx].set_color('r')
        self.pick = False

    def on_motion(self, event):
        """ Controls the sliding/stretching of the spectra """

        """
        Ignore this if we aren't in slide/stretch mode (ie pressing the
            mouse button
        """
        if self.domotion is False:
            return

        xdata = self.xdata.copy()
        """ Left mouse button is for stretching """
        if event.button==1 and (event.xdata is not None)  \
                and (event.ydata is not None) and event.xdata>self.lbound \
                and event.xdata<self.rbound:
            leftpts = self.xdata[self.lindx+1:self.indx+1].copy()
            left = scipy.linspace(leftpts[0],event.xdata,leftpts.size)
            rightpts = self.xdata[self.indx:self.rindx].copy()
            right = scipy.linspace(event.xdata,rightpts[-1],rightpts.size)[1:]
            xd = scipy.concatenate((left,right))
            xdata[self.lindx+1:self.rindx] = xd.copy()
            self.data.set_xdata(xdata)
        """ Middle mouse button is for sliding """
        if event.button==2:
            offset = event.xdata-self.origx
            xdata = xdata + offset
            self.data.set_xdata(xdata)
        for i in range(self.datapeaks.size):
            x = xdata[self.datapeaks[i]]
            l = self.datalines[i]
            l.set_xdata([x,x])
        plt.draw()

    def on_release(self, event):
        """ If the mouse button is released, reset! """
        if self.domotion:
            self.domotion = False
            self.xdata = self.spec.get_xdata().copy()
            plt.draw()

    def key_press(self,event):
        """
        m is for fitting, p is for selecting new lines (or leaving select
            mode) 
        """
        if type(event.key)!=type('m'):
            return
        if event.key.lower()=='m':
            self.do_fit()
        elif event.key.lower()=='d':
            self.do_fit(True)
        elif event.key.lower()=='a':
            if self.pick:
                self.pick = False
                self.canvas.mpl_disconnect(self.addid)
                self.connect()
                self.pick = False
            else:
                self.pick = True
                self.disconnect()
                print "Choose the line to add (a to exit)"
                self.addid = self.canvas.mpl_connect('button_press_event',
                                                        self.add_line)
        elif event.key.lower()=='w':
            self.write()
        elif event.key.lower()=='r':
            self.read()
        elif event.key.lower()=='n':
            self.data.remove()
            self.sky.remove()
            for i in self.lines:
                i.remove()
            a,b,c,d,e = self.start
            self.disconnect()
            self.canvas.mpl_disconnect(self.keyid)
            self.__init__(a,b,c,d,e)
            plt.draw()
        elif event.key.lower()=='q':
            if self.soln is None:
                self.do_fit()
            plt.close(self.data.get_figure())
            # stick a disconnect command here and see what happens...
            self.canvas.mpl_disconnect(self.keyid)

    def do_fit(self,fullFit=False):
        p = self.peaks[self.selected==1]
        fitdata = scipy.empty((self.xdata.size,2))
        fitdata[:,0] = scipy.arange(fitdata.shape[0])
        fitdata[:,1] = self.spec.get_xdata()
        fitdata = fitdata[(fitdata[:,0]>p[0]-20)&(fitdata[:,0]<p[-1]+20)]
        ord = int(raw_input('Enter order of fit: '))
        fit = sf.lsqfit(fitdata,'polynomial',ord)
        self.soln = fit
        if fullFit==True:
            from scipy import interpolate,optimize
            spec = self.spec.get_ydata()
            xvals = numpy.arange(spec.size).astype(numpy.float32)
            def opt(p):
                p = numpy.array(p)
                n = p[-1]
                coeff = numpy.atleast_2d(p[:-1]).T
                m = {'coeff':coeff,'type':'polynomial'}
                w = sf.genfunc(xvals,0.,m)
                mod = n*interpolate.splev(w,self.skymodel)
                return (spec-mod)/abs(spec)**0.5
            pars = fit['coeff'].flatten().tolist()
            pars.append(1.)
            coeff,ier = optimize.leastsq(opt,pars,maxfev=10000,epsfcn=1e-5)
            fit['coeff'] = numpy.atleast_2d(coeff[:-1]).T
            self.soln = fit
            print "pixel fit complete"
        else:
            print "Fit Complete"
        xdata = sf.genfunc(scipy.arange(self.xdata.size),0.,fit)
#        self.sky.set_xdata(xdata)
        self.data.set_xdata(xdata)
        for i in range(self.datapeaks.size):
            x = xdata[self.datapeaks[i]]
            l = self.datalines[i]
            l.set_xdata([x,x])
        plt.draw()
    def add_line(self,event):
        if self.canvas.toolbar.mode!='':
            if event.button==2:
                self.canvas.toolbar.zoom()
                self.canvas.toolbar.pan()
                self.canvas.toolbar.pan()
            return
        if event.xdata==None or event.ydata==None:
            print 'Invalid data'
            return
        xpos = event.xdata
        xdata = self.spec.get_xdata()
        ydata = self.spec.get_ydata()
        p = ndimage.maximum_filter(ydata,9)
        ind = abs(xdata-xpos).argmin()
        p = scipy.where((p==ydata))[0]
        if abs(p-ind).min()>5:
            print 'Not a line'
            return
        indx = p[abs(p-ind).argmin()]
        for i in self.peaks:
            if abs(indx-i)<9:
                print 'Too close to another line.'
                return
        if indx<5. or indx>xdata.size-6:
                print 'Too close to edge.'
                return
        peaks = scipy.arange(self.peaks.size+1)
        n = self.peaks[self.peaks<indx].size
        peaks[:n] = self.peaks[:n].copy()
        peaks[n] = indx
        peaks[n+1:] = self.peaks[n:].copy()
        sel = scipy.arange(peaks.size)
        sel[:n] = self.selected[:n].copy()
        sel[n] = 0
        sel[n+1:] = self.selected[n:].copy()
        self.peaks = peaks.copy()
        self.selected = sel.copy()
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        if self.spec==self.data:
            l = plt.axvline(xdata[indx],c='k')
            self.lines.insert(n,l)
            self.datapeaks = self.peaks
            self.datasel = self.selected
            self.datalines = self.lines
        else:
            l = plt.axvline(xdata[indx],c='k',ls=':')
            self.lines.insert(n,l)
            self.skypeaks = self.peaks
            self.skysel = self.selected
            self.skylines = self.lines
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        plt.draw()
        self.canvas.mpl_disconnect(self.addid)
        self.connect()
        self.pick = False
    def disconnect(self):
        self.canvas.mpl_disconnect(self.pressid)
        self.canvas.mpl_disconnect(self.moveid)
        self.canvas.mpl_disconnect(self.offid)
    def write(self):
        oname = raw_input('Name of output file: ')
        f = open(oname,'w')
        xdata = self.spec.get_xdata()
        selected = self.datasel
        peaks = self.datapeaks
        soln = self.soln
        pickle.dump([xdata,selected,peaks,soln],f)
        f.close()
        print "Writing Complete"
    def read(self):
        oname = raw_input('Name of input file: ')
        try:
            f = open(oname,'r')
            xdata,selected,peaks,soln = pickle.load(f)
            f.close()
        except:
            print "Could not read file: %s"%oname
            return
        self.data.set_xdata(xdata)
        self.spec.set_xdata(xdata)
        self.datapeaks = peaks
        self.datasel = selected
        self.soln = soln
        for i in self.lines:
            i.remove()
        for i in range(len(self.lines)):
            del self.lines[0]
        for i in range(len(peaks)):
            c = 'k'
            if selected[i]==1:
                c = 'r'
            l = plt.axvline(xdata[peaks[i]],c=c)
            self.datalines.append(l)
        self.lines = self.datalines
        self.peaks = self.datapeaks
        self.selected = self.datasel
        plt.draw()



def id_spec(spec,model):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fmt_xdata = plt.FormatStrFormatter('%4.2f')
    ax.fmt_ydata = plt.FormatStrFormatter('%4.2f')
    from scipy import ndimage,interpolate

    skymodel = model['model']
    data = spec.copy()
    blue,red,scale = model['blue'],model['red'],model['scale']
    while (red-blue)/scale<data.size:
        red += scale
        blue -= scale
    wave = scipy.arange(blue,red,scale)
    blue,red = model['blue'],model['red']
    sky = interpolate.splev(wave,skymodel)

    sky /= sky.mean()/data.mean()
    plt.plot(wave,sky,c='gray')

    from scipy import optimize,signal
    corr = signal.correlate(sky,data,mode='valid')
    w0 = wave[corr.argmax()]
    p = [w0,scale,0.,0.,1.]
    xvals = numpy.arange(data.size).astype(numpy.float32)
    def opt(p):
        p = numpy.array(p)
        n = p[-1]
        coeff = numpy.atleast_2d(p[:-1]).T
        m = {'coeff':coeff,'type':'polynomial'}
        w = sf.genfunc(xvals,0.,m)
        mod = n*interpolate.splev(w,skymodel)
        return (data-mod)/abs(data)**0.5
    #coeff,ier = optimize.leastsq(opt,p,maxfev=10000,epsfcn=1e-5)
    #fit = {'coeff':numpy.atleast_2d(coeff[:-1]).T,'type':'polynomial'}
    #dwave = sf.genfunc(xvals,0.,fit)
    dwave = numpy.linspace(blue,red,data.size)
    dr = IDSpectrum(data,sky,dwave,wave,skymodel)
    print "dr.keyid:", dr.keyid
    print "dr.soln:", dr.soln

    x = numpy.arange(sky.size)
    w = sf.genfunc(x,0.,dr.soln)

    if dr.soln is None:
        dr.do_fit()
    xdata = dr.data.get_xdata()
    ydata = dr.data.get_ydata()
    return dr.soln
