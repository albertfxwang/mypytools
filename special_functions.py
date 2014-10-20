import scipy
from scipy import optimize

#
# Basic (polynomial) fitting routine frontends
#
def lsqfit(input,fittype,xorder,yorder=0):
    p = scipy.ones((xorder+1,yorder+1))
    par = {'coeff':p,'type':fittype}
    return lsqfitter(input,par)

def lsqfitter(input,par):
    par = par.copy()
    p = par['coeff']
    # If only one array is input, assume they are z-values
    if input.ndim==1:
        input = scipy.atleast_2d(input).T
    xorder = p.shape[0]-1
    yorder = p.shape[1]-1
    if input.shape[1]>3:
        print "Incorrect input array size: should be nx2 or nx3!"
        return 0
    elif input.shape[1]==1:
        x = scipy.arange(0,input.shape[0])
        y = 0
        z = input[:,0].copy()
    elif input.shape[1]==2:
        x = input[:,0]
        y = 0
        z = input[:,1]
    else:
        x = input[:,0]
        y = input[:,1]
        z = input[:,2]

    good = scipy.isfinite(z)
    x = x[good]
    z = z[good]

    try:
        y = y[good]
    except:
        pass

    coeffs = unpack_coeff(par)

    fit,ier = optimize.leastsq(dolsqfit,coeffs,args=(x,y,z,par),maxfev=6000)

    return build_coeff(fit,par)

# Convert from leastsq par array to coeffs array
def build_coeff(p,par):
    xorder = par['coeff'].shape[0]-1
    yorder = par['coeff'].shape[1]-1
    pars = scipy.zeros(par['coeff'].shape)
    if xorder<yorder:
        order = xorder
    else:
        order = yorder
    ncoeffs = xorder+1+yorder+((order*(order-1))/2)
    k = 0
    for i in range(yorder+1):
        for j in range(xorder+1):
            if i+j>order and i>0 and j>0:
                break
            pars[j,i] = p[k]
            k += 1
    out = {'coeff':pars,'type':par['type']}
    return out

# Convert from coeffs array to leastsq par array
def unpack_coeff(par):
    p = par['coeff'].copy()
    xorder = p.shape[0]-1
    yorder = p.shape[1]-1

    if xorder<yorder:
        order = xorder
    else:
        order = yorder
    ncoeffs = xorder+1+yorder+((order*(order-1))/2)
    coeffs = scipy.zeros(ncoeffs)

    k = 0
    for i in range(yorder+1):
        for j in range(xorder+1):
            if i+j>order and i>0 and j>0:
                break
            coeffs[k] = p[j,i]
            k += 1
    return coeffs


# Fitting function for leastsq
def dolsqfit(coeffs,x,y,z,par):
    p = build_coeff(coeffs,par)
    return z - genfunc(x,y,p)


# Generic curve/surface making routine (returns an array)
#   This is about 2x faster than using the scipy.special package.
def genfunc(x,y,par):
    p = scipy.atleast_2d(par['coeff'])
    type = par['type']
    x = scipy.atleast_1d(x)
    y = scipy.atleast_1d(y)

    # If fitting a constant x or y to a two-d fit, make the correct array
    if p.ndim==2:
        if x.size==1:
            xtmp = x[0]
            x = y*0. + xtmp
        elif y.size==1:
            ytmp = y[0]
            y = x*0. + ytmp

    value = scipy.zeros(x.size)

    # These arrays can end up being huge, so only do NPOINTS points
    #   at a time. After a lot of time trials, segments of 5000 are about
    #   the optimal size (a factor of 2 faster than 2000000, for example).
    NPOINTS = 5000
    start = 0
    while start<x.size:
        end = start+NPOINTS
        if end>x.size:
            end = x.size
        len = end-start
        xfit = x[start:end]
        # Ensure we are fitting in both dimensions!
        if y.size==x.size:
            yfit = y[start:end]
        p_x = scipy.zeros((p.shape[0],len))
        p_y = scipy.zeros((p.shape[1],len))
        for i in range(p.shape[0]):
            if i==0:
                p_x[i] = 1.
            elif i==1:
                p_x[i] = xfit
                if type=='hermite':
                    p_x[i] *= 2.
            else:
                if type=="legendre":
                    a = (2.*i-1.)/i
                    c = 1./i
                elif type=="hermite":
                    a = 2.
                    c = 2.*i-2.
                elif type=="chebyshev":
                    a = 2.
                    c = 1.
                elif type=="polynomial":
                    a = 1.
                    c = 0.
                else:
                    a = 1.
                    c = 0.
                p_x[i] = a*xfit*p_x[i-1]-c*p_x[i-2]

        for i in range(p.shape[1]):
            if i==0:
                p_y[i] = 1.
            elif i==1:
                p_y[i] = yfit
                if type=='hermite':
                    p_y[i] *= 2.
            else:
                if type=="legendre":
                    a = (2.*i-1.)/i
                    c = 1./i
                elif type=="hermite":
                    a = 2.
                    c = 2.*i-2.
                elif type=="chebyshev":
                    a = 2.
                    c = 1.
                elif type=="polynomial":
                    a = 1.
                    c = 0.
                else:
                    a = 1.
                    c = 0.
                p_y[i] = a*yfit*p_y[i-1] - c*p_y[i-2]

        if p.shape[1]>p.shape[0]:
            order = p.shape[0]-1
        else:
            order = p.shape[1]-1
        for i in  range(p.shape[1]):
            for j in range(p.shape[0]):
                if i+j>order and i>0 and j>0:
                    break
                value[start:end] += p[j,i]*p_x[j]*p_y[i]
        start = end

    return value


# Fit "n" 1d gaussians!
def ngaussfit(data,p,weight=0):
    return nmodelfit(data,p,"gauss",weight)


def nmodelfit(data,p,model,weight=0):
    if p.ndim==2:
        if p.shape[0]==2:
            p = p.T
        mask = p[:,1].copy()
        t = scipy.zeros(mask.sum())
        static = scipy.zeros(mask.size-t.size)
        j = 0
        k = 0
        for i in range(mask.size):
            if mask[i]>0:
                t[j] = p[i,0]
                j += 1
            else:
                static[k] = p[i,0]
                k += 1
        p = t.copy()
    else:
        mask = scipy.ones(p.size)
        static = scipy.zeros(0)

    if data.ndim==1:
        x = scipy.arange(0.,data.size,1.)
        z = data.copy()
    else:
        x = data[:,0]
        z = data[:,1]

    good = scipy.isfinite(z)
    x = x[good]
    z = z[good]

    pars,cov,info,mesg,ier = optimize.leastsq(domodel,p,args=(x,z,mask,static,model,weight),maxfev=10000,full_output=True)

    chi2 = info['fvec']
    chi2 = chi2*chi2
    if weight==0:
        chi2 /= abs(z)
    chi2 = chi2.sum()

    p = scipy.zeros(mask.size)
    j = 0
    k = 0
    for i in range(mask.size):
        if mask[i]>0:
            p[i] = pars[j]
            j += 1
        else:
            p[i] = static[k]
            k += 1
        if i%3==0:
            p[i] = scipy.fabs(p[i])
    return p,chi2

def dogauss(p,x,z,mask,static):
    par = scipy.zeros(mask.size)
    j = 0
    k = 0
    for i in range(mask.size):
        if mask[i]>0:
            par[i] = p[j]
            j += 1
        else:
            par[i] = static[k]
            k += 1
    model = ngauss(x,par)
    diff = z - model
    return diff

def domodel(p,x,z,mask,static,model,weight):
    par = scipy.zeros(mask.size)
    j = 0
    k = 0
    for i in range(mask.size):
        if mask[i]>0:
            par[i] = p[j]
            j += 1
        else:
            par[i] = static[k]
            k += 1
    if model=="gauss":
        vals = ngauss(x,par)
    elif model=="moffat":
        vals = nmoffat(x,par)
    sigma = scipy.sqrt(abs(vals))
    sigma[sigma<1e-7] = 1.e-7
    if weight==0:
        sigma = 1.
    diff = (z-vals)/sigma
    return diff


def ngauss(x,p):
    from math import fabs
    n = (p.size-1)/3 # Number of gaussians being fit
    value = scipy.ones(x.size)
    value *= p[0]
    for i in range(n):
        k = i*3
        eval = (x-p[k+2])/p[k+3]
        eval = eval*eval/-2.
        value += p[k+1]*scipy.exp(eval)
        p[k+3] = fabs(p[k+3])

    return value


def nmoffat(x,p):
    n = (p.size-1)/4
    value = scipy.zeros(x.size)
    value += p[0]
    for i in range(n):
        k = i*4
        amp = p[k+1]
        cent = p[k+2]
        fwhm = p[k+3]
        b = scipy.fabs(p[k+4])
        a = 0.5*fwhm/scipy.sqrt(scipy.power(2.,1./b)-1.)
        r = x-cent
        value += amp*scipy.power(1.+(r/a)*(r/a),-1.*b)
        p[k+3] = scipy.fabs(p[k+3])
    return value
