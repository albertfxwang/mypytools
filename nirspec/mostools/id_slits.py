import scipy,numpy
from scipy import ndimage,signal,interpolate


def clip(arr,thresh=3.5):
    """
    clip(arr,thresh=3.5)

    Simple sigma-clipping algorithm. Returns avg,std of clipped array.
    """
    a = numpy.array(arr)

    avg,std = a.mean(),a.std()
    while 1:
        avg,std,size = a.mean(),a.std(),a.size
        a = a[abs(a-avg)<thresh*std]
        if size==a.size:
            break
    return avg,std


def Clip(arr,nsig=3.5,hicut=0.9,locut=0.1):
    a = numpy.array(arr)
    a.sort()
    l = a.size
    a = a[locut*l:hicut*l]
    m,s,l = a.mean(),a.std(),a.size
    while 1:
        a = a[abs(a-m)<nsig*s]
        if a.size==l:
            return m,s
        m,s,l = a.mean(),a.std(),a.size


def check_star(peaks,data):
    """
    check_star(peaks,data)

    Determines whether or not a slit looks like it is a starbox. This is
      done by simply checking the 3pixels bordering each peak and ensuring
      that none are less than half of the peak (ie that the FWHM>7 pixels).

    Returns True if more than half of the peaks look like boxes, otherwise
      returns False.
    """
    star = 0
    for i in peaks:
        max = data[i]
        if i<3 or i+4>data.size:
            continue
        mean = data[i-3:i+4].mean()
        if (max-mean)<0.1*max:
            star += 1
    if star*2>peaks.size:
        return True
    else:
        return False


def findlines(z,bgsub=True,SATURATED=57000.):
    """
    findlines(z)

    Quickly find the peaks of arclines. Returns a list containing the peak
      locations.
    """
    z = z.copy()
    s = z.copy()

    """ First identify peaks. """
    max = ndimage.maximum_filter(z,9)
    p = scipy.where((max==z)&(z<SATURATED)&(max>0))[0]
    s = z[p]

    """ Reject low peaks. """
    bg = ndimage.percentile_filter(s,10,21)
    peaks = scipy.where(s>bg*5.)[0]
    return p[peaks]


def id_slits(arc,find_stars=True,chilimit=2.5,SATURATED=57000.,useLines=True):
    """
    id_slits(arc,find_stars=True)

    Determine the top/bottom of each slit in the 2d y-corrected arc image.

    If find_stars is True, returns slit,starbox, otherwise returns slit.
    """

    arc = arc.copy()
    """ Attempt to avoid saturated lines """
    w = arc.shape[1]
    tmp = arc.copy()
    tmp[tmp>SATURATED] = 0.
    tmpSorted = scipy.sort(tmp,axis=1)
    flux = tmpSorted[:,w*0.97:w*0.98].mean(axis=1)
    minflux = scipy.median(flux)/4.
    del tmp

    if find_stars:
        starbox = []
    slit = []

    if useLines==False:
        flux = scipy.sort(arc,1)[:,w*4/5]
        minflux = scipy.median(flux[flux.size/3:flux.size*2/3])/2.
        mask = scipy.where(flux>minflux,1.,0.)
        inSlit = False
        tmp = []
        meds = []
        for i in range(mask.size):
            if inSlit:
                if mask[i]==0:
                    inSlit = False
                    end = i-1
                    if end-start>8:
                        tmp.append([start+1,end-1])
                        slit = arc[start+3:end-3,100:-100].mean(0)
                        meds.append(slit.max())
            elif mask[i]==1:
                start = i
                inSlit = True
        if inSlit:
            end = i
            if end-start>8:
                tmp.append([start+1,end-1])
                slit = arc[start+3:end-3,100:-100].mean(0)
                meds.append(slit.max())
        meds = numpy.array(meds)
        if find_stars:
            slit = []
            starbox = []
            m,s = Clip(meds,nsig=3.,locut=0.,hicut=0.75)
            for i in range(len(tmp)):
                if meds[i]<m+s*5:
                    slit.append(tmp[i])
                else:
                    starbox.append(tmp[i])
            return slit,starbox
        return tmp

        m,s = clip(tmpSorted[arc.shape[0]/2,:w*0.05],2.)

    inSlit = False
    i = 0
    while i<arc.shape[0]:
#        lines = findlines(arc[i])
        if useLines:
            lines = findlines(arc[i])
        else:
            med = scipy.median(arc[i])
            if med>m+5*s:
                lines = [0]*10
            else:
                lines = [0]
        if len(lines)<9 and inSlit==False:
            i += 1
            continue
        elif len(lines)>9 and inSlit==False:
            inSlit = True
            start = i
            i += 1
            continue

        bestchi = 1e29
        if len(lines)>9:
            #bestchi = 1e29
            x = scipy.arange(arc[i].size)
            smooth = ndimage.gaussian_filter(arc[i],1.)
            model = interpolate.splrep(x,smooth,k=3,s=0)
            comp = ndimage.gaussian_filter(arc[i-1],1.)
            usedpix = arc[i-1][10:-10]>scipy.median(arc[i-1])
            for o in range(30):
                offset = float(o-15.)/5.
                row = interpolate.splev(x[10:-10]+offset,model)
                chi = (comp[10:-10]-row)**2/(abs(comp[10:-10]))
                chi = chi[usedpix]
                chi.sort()
                chi = chi[:-chi.size/100] # Reject the five highest points
                if chilimit>6. and i>600 and o>6 and 1==2:
                        import pylab
                        pylab.plot(row)
                        pylab.plot(comp[10:-10])
                        pylab.figure()
                        pylab.plot((row-comp[10:-10])**2/(abs(comp[10:-10])+16.))
                        pylab.show()
                if chi.sum()/chi.size<bestchi:
                    bestchi = chi.sum()/chi.size

        if inSlit is True and (bestchi>chilimit or len(lines)<9):
            """ The row is at the top edge of the slit. """
            inSlit = False
            end = i

            i += 1
            if end-start<3:
                continue
            """
            Conservatively shrink the edges. A better approach
              might be to use the flatfield data and set the edge
              to where the flux is, say, 1 sigma below the nominal
              level for the slit.
            """
#            if start!=0:
#                start += 2
#            end -= 2

            """ Check if the slit is a starbox (if requested) """
            if find_stars:
                mid = (start+end)/2
                peaks = findlines(arc[mid],False)
                is_star = check_star(peaks,arc[mid])
            else: is_star = False

            """ Skip small slits """
            if not is_star and end-start<11:
                continue
            elif is_star and end-start<9:
                continue

            """
            Conservatively shrink the edges. A better approach
              might be to use the flatfield data and set the edge
              to where the flux is, say, 1 sigma below the nominal
              level for the slit.
            """
            if is_star:
                starbox.append([start,end])
            else:
                while flux[start+1]-flux[start]>3.*flux[start]**0.5:
                    start += 1
                while flux[end-1]-flux[end]>3.*flux[end]**0.5:
                    end -= 1
                if flux[start:end].mean()<minflux:
                    continue
                slit.append([start,end])

        elif i+1==arc.shape[0] and end<start:
            """ The top of the mask is also the top of a slit. """
            end = i+1
            if find_stars:
                mid = (start+end)/2
                peaks = findlines(arc[mid],False)
                is_star = check_star(peaks,arc[mid])
            else: is_star = False

            if not is_star and end-start<11:
                continue
            elif is_star and end-start<9:
                continue

            if is_star:
                starbox.append([start+2,end])
            else:
                while flux[start+1]-flux[start]>3.*flux[start]**0.5:
                    start += 1
                if flux[start:end].mean()<minflux:
                    continue
                slit.append([start,end])
            break
        else:
            """ In the middle of the slit, nothing to do.... """
            i += 1

    if find_stars:
        return slit,starbox
    return slit
