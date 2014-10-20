# Cosmology tools

from numpy import *

C=2.99792458e10
SEC_PER_YEAR=3.15569e7
CM_PER_MPC=3.08568025e24
CM3_PER_MPC3=(CM_PER_MPC*CM_PER_MPC*CM_PER_MPC)
SKYARCMIN       = 1.4850936e8   # arcmin in 4pi steradian


# Global variables
iv_H0 = 0.  # For integrating volume
iv_q0 = 0.  # For integrating volume
sl_H0 = 0. 
sl_q0 = 0.
sl_omega_l = 0.
sl_omega_k = 0.
sl_omega_m = 0.

def h(z,omega_tot,omega_m,omega_l):
    return sqrt(omega_l+(1-omega_l-omega_tot)*(1+z)**2 + omega_tot*(1+z)**3)


def H0_hertz(H0):
    return H0*3.240777649e-20

def lumdist(z,H0,omega_tot,omega_m,omega_l):
    """ Calculate luminosity distance for a flat universe"""
    return (1+z)*pmdist(z,H0,omega_tot,omega_m,omega_l)

def pmdist(z,H0,omega_tot,omega_m,omega_l):
    if type(z) == type(arange(10.)):
        pmd = 0.*z
        for i in range(len(z)):
            pmd[i] = pm_distance(z[i],H0,omega_tot,omega_m,omega_l)
        return pmd
    else:
        return pm_distance(z,H0,omega_tot,omega_m,omega_l)

def luminosity_distance(z,H0,omega_tot,omega_m,omega_l):
    dm= pm_distance(z,H0,omega_tot,omega_m,omega_l)
    ld = (1+z)*dm
    return ld


def pm_distance(z,H0,omega_tot,omega_m,omega_l):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0

    H0 = H0_hertz(H0)
    sl_H0   = H0 
    sl_omega_m = omega_m
    sl_omega_k = 1. - omega_tot
    sl_omega_l = omega_l

    if z < 0.0001:
        dm = z * C / ((1+z)*H0)        # cm/s 
        return dm

    ok = sqrt(abs(sl_omega_k))
#   lzsamp = arange(-40.,log10(z),(log10(z)+40)/10000.)
#   zsamp = arange(0.,z,z/1000.)
#   zsamp = 10.**lzsamp
#   e = aa(zsamp)
#   integral = trapz(e,x=zsamp)
    integral = closedpoints(aa,0,z,TOL=1.e-9)

    if  sl_omega_k > 0:
            dm = (C/(ok*sl_H0))*sinh(ok*integral)
    if  sl_omega_k < 0:
            dm = (C/(ok*sl_H0))*sin(ok*integral)
    if  sl_omega_k == 0:
            dm = (C/sl_H0)*integral
    return dm

def da12(z1,z2,H0,omega_tot,omega_m,omega_l):
    "da12(z1,z2,H0,omega_tot,omega_m,omega_l) - returns angular-diameter distance between objects at z1,z2 in cm"
    omega_k = 1 - (omega_m+omega_l)
    dm1 = pmdist(z1,H0,omega_tot,omega_m,omega_l)
    dm2 = pmdist(z2,H0,omega_tot,omega_m,omega_l)
    dh = C/H0_hertz(H0)
    term1 = (1.0/(1.0+z2))
    term2 = dm2*sqrt(1.0+omega_k*dm1**2/dh**2)
    term3 = dm1*sqrt(1.0+omega_k*dm2**2/dh**2)
    da = term1*(term2-term3)
    return(da)

def cm_to_mpc(d):
    return d/CM_PER_MPC

def cm_to_distmod(d):
    """distmod(d) - returns distance modulus from distance in cm"""
    return(5.0*log10(cm_to_mpc(d))+25.0)

def distmod(z,H0,omega_tot,omega_m,omega_l):
    """distmod(d) - returns distance modulus corresponding to luminosity-dist"""
    return (cm_to_distmod(lumdist(z,H0,omega_tot,omega_m,omega_l)))

def lookback(z,H0,omega_tot,omega_m,omega_l):
    if type(z) == type(arange(10.)):
        lbt = 0.*z
        for i in range(len(z)):
            lbt[i] = lookbacktime(z[i],H0,omega_tot,omega_m,omega_l)
        return lbt
    else:
        return lookbacktime(z,H0,omega_tot,omega_m,omega_l)

def lookbacktime(z,H0,omega_tot,omega_m,omega_l):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    sl_H0   = H0_hertz(H0)
    sl_omega_m = omega_m
    sl_omega_k = 1 - omega_tot
    sl_omega_l = omega_l
    lb = closedpoints(dt_z,0.,z)
    lb_time = lb / SEC_PER_YEAR
    return lb_time

def age_from_z(z,zform,H0,omega_tot,omega_m,omega_l):
   t = lookback(z,H0,omega_tot,omega_m,omega_l)
   tf = lookback(zform,H0,omega_tot,omega_m,omega_l)
   tg = tf - t
   return tg

def volume_element(z,H0,omega_tot,omega_m,omega_l):
    z1 = z*1.1;
    z2 = z/1.1;
    v1 = integrate_volume(z1,H0,omega_tot,omega_m,omega_l)
    v2 = integrate_volume(z2,H0,omega_tot,omega_m,omega_l)
    dvdz = (v2-v1)/(z2-z1);
    return dvdz

def iv_volume(dm):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    dv = dm*dm/sqrt(1+sl_omega_k*sl_H0*sl_H0*dm*dm/(C*C))
    dv *= 4*pi / CM3_PER_MPC3
    return dv 

def integrate_volume(z,H0,omega_tot,omega_m,omega_l):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    dm = pm_distance(z,H0,omega_tot,omega_m,omega_l)
    sl_H0 = H0_hertz(H0)
    sl_omega_k = 1-omega_tot;
    v = closedpoints(iv_volume,0.,dm,TOL=1.e-9);
    return v

def dt_z(z):
    """ Compute time element """
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    x1 = aa(z)
    x2 = 1./(sl_H0 * (1+z))
    dtdz = x1*x2
    return dtdz

def aaa(z):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    x = 1+z
    a = x/sqrt(x*x*(1+sl_omega_m*z) - z*(2+z)*sl_omega_l)
    return a 

def aa(z):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    x = 1+z
    a = 1./sqrt(x*x*(1+sl_omega_m*z) - z*(2+z)*sl_omega_l)
    return a

def E(z):
    global sl_omega_l
    global sl_omega_k
    global sl_omega_m
    global sl_H0
    x = 1+z
    h2 = (sl_omega_m*x*x*x+1-sl_omega_m)
    return 1/sqrt(h2)
    
def dv(z1,z2,H0,omega_tot,omega_m,omega_l):
    volume1 = volume(z1,H0,omega_tot,omega_m,omega_l)
    volume2 = volume(z2,H0,omega_tot,omega_m,omega_l)
    return volume2 - volume1

def dv_arcmin(z1,z2,H0,omega_tot,omega_m,omega_l):
    v = dv(z1,z2,H0,omega_tot,omega_m,omega_l)
    return(v/SKYARCMIN)

def volume_arcmin(z,H0,omega_tot,omega_m,omega_l):
    y = volume(z,H0,omega_tot,omega_m,omega_l)
    return(y/SKYARCMIN)

def volume(z,H0,omega_tot,omega_m,omega_l):
    if type(z) == type(arange(10.)):
        vol = 0.*z
        for i in range(len(z)):
            vol[i] = integrate_volume(z[i],H0,omega_tot,omega_m,omega_l)
        return vol 
    else:
        return  integrate_volume(z,H0,omega_tot,omega_m,omega_l)


def closedpoints(func, a, b, TOL=1e-6):		# f(x)=func(x)
    """
    Closed Simpson's rule for 
        \int_a^b f(x) dx
    Divide [a,b] iteratively into h, h/2, h/4, h/8, ... step sizes; and,
    for each step size, evaluate f(x) at a+h, a+3h, a+5h, a+7h, ..., b-3h,
    b-h, noting that other points have already been sampled.
    
    At each iteration step, data are sampled only where necessary so that
    the total data is represented by adding sampled points from all
    previous steps:
        step 1:	h	a---------------b
        step 2:	h/2 	a-------^-------b
        step 3:	h/4	a---^-------^---b
        step 4:	h/8	a-^---^---^---^-b
        total:		a-^-^-^-^-^-^-^-b
    So, for step size of h/n, there are n intervals, and the data are
    sampled at the boundaries including the 2 end points.
    
    If old = Trapezoid formula for an old step size 2h, then Trapezoid
    formula for the new step size h is obtained by 
        new = old/2 + h{f(a+h) + f(a+3h) + f(a+5h) + f(a+7h) +...+ f(b-3h)
	    + f(b-h)}
    Also, Simpson formula for the new step size h is given by
        simpson = (4 new - old)/3
    """
    h = b - a
    old2 = old = h * (func(a) + func(b)) / 2.0
    count = 0
    while 1:
	h = h / 2.0
	x, sum = a + h, 0
	while x < b:
	    sum = sum + func(x)
	    x = x + 2 * h
	new = old / 2.0 + h * sum
	new2 = (4 * new - old) / 3.0
	if abs(new2 - old2) < TOL * (1 + abs(old2)): return new2
	old = new	# Trapezoid
	old2 = new2	# Simpson
	count = count + 1
#	print 'closedpoints(%d): Trapezoid=%s, Simpson=%s' % (count, new, new2)
