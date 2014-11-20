#!/usr/bin/python
from numpy import *
import os, popen2, subprocess

### contains frequently used tasks using BC03 code

### define a function that generates composite stellar population from BC03 code
def gencsp(sspname,tau,cspname,dust='n',sf_hist=1,rcy='n',rcy_frac=0,max_age=20):
    '''MUST OPEN THE PIPELINE TO CSP_GALAXEV FIRST!!!
       sspname: input ssp file name
       tau: star formation time scale (in Gyr)
       cspname: output csp file name
       dust: include dust extinction or not, 'n' is no, 'y' is yes
       sf_hist: 1 = exponential
       rcy: recyle gas or not
       rcy_frac: recycled fraction if rcy = 'y'
       max_age: maximum time after which flux = 0'''
    #if not csp: raise 'pipeline not open yet'
    csp = os.popen('/devel/goods14/khuang/bc03/src/csp_galaxev','w')
    print >> csp, sspname
    print >> csp, dust    #include dust or not
    print >> csp, sf_hist   #decide star formation history
    print >> csp, tau     #star formation time scale
    print >> csp, rcy     #recyle gas or not
    if rcy == 'y':
        print >> csp, rcy_frac
    print >> csp, max_age  #maximum age
    print >> csp, cspname
    csp.close()

def genmod(cspname, age, sedname, wavemin=91.0, wavemax=1.6e6):
    #sout, gpl = popen2.popen2('/devel/goods14/khuang/bc03/src/galaxevpl')
    f = subprocess.Popen(['/devel/goods14/khuang/bc03/src/galaxevpl'],stdin=subprocess.PIPE)
    print >> f.stdin, "%s" % cspname
    print >> f.stdin, '%f, %f' % (wavemin, wavemax)   # enter the desired wavelength range
    print >> f.stdin, "%f" % age   # should be in Gyr
    print >> f.stdin, "%s" % sedname
    #retcode = f.wait()
    #print retcode

def read4color(fname):
    if fname[-6:] != '4color':
        raise NameError, 'input file not a *.4color file'
    f = open(fname, 'r')
    ncolumn = 10
    #nrow = 0
    dlist = []   # dlist will hold all the data
    for line in f.readlines():
        if line[0] == '#':
            pass
        else:
            g = array(line.split())
            g = g.astype(float64)
            g = reshape(g, (ncolumn, 1))  # reshape for concatenating 
            if len(dlist) == 0:
                dlist = g   # create the first row of dlist
            else:
                dlist = concatenate((dlist,g),1)
    f.close()
    return dlist

def read_stmass(fname, logage):
    if fname[-6:] != '4color':
        raise NameError, 'input file not a *.4color file'
    dlist = read4color(fname)
    logage_array = dlist[0]
    stmfrac = dlist[6]   # stmfrac = stellar mass fraction of 1 M_solar
    # now interpolate the stellar mass
    age_index = searchsorted(logage_array, logage)
    logage_top = logage_array[age_index]
    stmfrac_top = stmfrac[age_index]
    logage_bottom = logage_array[age_index-1]
    stmfrac_bottom = stmfrac[age_index-1]
    age_top = 10**logage_top
    age_bottom = 10**logage_bottom   # convert from logage to age
    age = 10**logage
    intstmfrac = stmfrac_bottom + (stmfrac_top - stmfrac_bottom)*(age - age_bottom)/(age_top - age_bottom) 
    return intstmfrac 

def read_sfr(fname, logage):
    if fname[-6:] != '4color':
        raise NameError, 'input file not a *.4color file'
    dlist = read4color(fname)
    logage_array = dlist[0]
    sfrfrac = dlist[9]   # sfrfrac = star formation rate per total 1 M_solar
    # now interpolate the stellar mass
    age_index = searchsorted(logage_array, logage)
    logage_top = logage_array[age_index]
    sfrfrac_top = sfrfrac[age_index]
    logage_bottom = logage_array[age_index-1]
    sfrfrac_bottom = sfrfrac[age_index-1]
    age_top = 10**logage_top
    age_bottom = 10**logage_bottom   # convert from logage to age
    age = 10**logage
    intsfrfrac = sfrfrac_bottom + (sfrfrac_top - sfrfrac_bottom)*(age - age_bottom)/(age_top - age_bottom) 
    return intsfrfrac 

def read_mtolratio_b(fname, logage):
    if fname[-6:] != '4color':
        raise NameError, 'input file not a *.4color file'
    dlist = read4color(fname)
    logage_array = dlist[0]
    mtol_b = dlist[4]   # sfrfrac = star formation rate per total 1 M_solar
    # now interpolate the stellar mass
    age_index = searchsorted(logage_array, logage)
    logage_top = logage_array[age_index]
    mtol_b_top = mtol_b[age_index]
    logage_bottom = logage_array[age_index-1]
    mtol_b_bottom = mtol_b[age_index-1]
    age_top = 10**logage_top
    age_bottom = 10**logage_bottom   # convert from logage to age
    age = 10**logage
    intmtol_b = mtol_b_bottom + (mtol_b_top - mtol_b_bottom)*(age - age_bottom)/(age_top - age_bottom) 
    return intmtol_b

