#!/usr/bin/python
from numpy import *
import numpy as np
import os, popen2, subprocess, glob
import pysynphot as S
from scipy import integrate, interpolate
from stats import gauss
from pygoods import sextractor

L_solar = 3.826e33   # erg/s
pc_cm = 3.0857e18  # cm
area_10pc = 4 * np.pi * (10 * pc_cm)**2   # in cm^2
normfactor = L_solar / area_10pc  # a scaling factor from BC03 flux to erg/s/cm^2/A
h = 6.636e-27  # Planck constant; erg * s
# LyC_H is the conversion factor from the LyC photon counts to Hydrogen 
# recombination line strengths (in erg/s). Each line is specified by 
# its wavelength in Angstrom. The coefficients are taken from Leitherer &
# Heckman 1995, equations 9-12.
# NOTE: the first line is always H-beta, because we need it to scale the
# metal lines.
LyC_H = [[4861., 4.76e-13], [6563., 1.36e-12], [1282., 7.13e-14], 
         [2165., 1.31e-14]]
# The nebular line ratios relative to H-beta, given by Anders & 
# Fritze-v. Alvensleben 2003. There are different values for different 
# metallicities: Z1=0.002 Z_solar, Z2=0.02 Z_solar, Z3=0.4 Z_solar, 
# Z4=Z_solar, Z5=2.5 Z_solar. Values for Z >= Z3 are the same.
nebular_ratios = {'lambda':[1335., 1663., 1909., 2141., 2326., 2798., 3727.,
                           3869., 3889., 3970., 4026., 4068.6, 4076.35,
                           4363., 4471., 4711., 4958.91, 5006.84, 5199.,
                           5755., 5876., 6300., 6312., 6548.05, 6583.45,
                           6678.0, 6716.0, 6730.0, 7065.0, 7135.79, 7319.99,
                           7330.73, 7751.11, 9068.6, 9530.85, 10286.73,
                           10320.49, 10336.41],
                  'Z3':    [0.110, 0.010, 0.180, 0.010, 0.290, 0.070, 3.010,
                           0.300, 0.107, 0.159, 0.015, 0.029, 0.011,
                           0.010, 0.050, 0.0, 1.399, 4.081, 0.030,
                           0.010, 0.140, 0.130, 0.030, 0.136, 0.404,
                           0.030, 0.300, 0.210, 0.040, 0.035, 0.026,
                           0.014, 0.086, 0.945, 0.365, 0.048,
                           0.058, 0.054]}

# Now also implement nebular continuum emission, mostly from HI, HeII, HeI, 
# and the 2-photon process. The coefficients are taken from Brown & Mathews
# 1970, which I followed from Krueger et al. 1995, and which was referenced 
# by Anders et al. 2003.
# The Flam component from the nebular continuum emission was from equation (7)
# and (8) of Krueger et al. 1995.
# At 10000 K; in units of 1e-40 erg/cm^2/s/Hz; from Table 1 of BM1970
gamma_HI = [[2500.,2600.,3122.,3646.,4000.,4500.,5696.,7000.,8204.,10000.,10100.],
            [0.,5.39,13.19,1.380,1.946,2.88,5.52,8.67,4.31,5.86,0.]]
n_HeII = 0.  # from Krueger et al. 1995; so ignore HeII at 10000 K...
gamma_HeI = [[2500.,2600.,3122.,3422.,3680.,4000.,4500.,5696.,6636.,7440.,7849.,8195.,8197.,8268.,10000.,10100.],
            [0.,11.87,17.45,23.0,5.69,2.07,3.31,6.76,9.09,10.18,10.61,8.50,5.88,5.10,5.86,0.]]
n_HeI = 0.095  # from Krueger et al. 1995
# The two-photon decay emission; Table 4 of Brown & Mathews 1970
gamma_2q = [[2300.,2431.,2701.,3039.,3473.,4052.,4863.,6078.,8104.,12157.,24313.,25000.],
            [0.,7.22,6.46,5.68,4.84,4.00,3.16,2.33,1.536,0.818,0.254,0.]]

def calc_Lyc(specfile, unit='Lsolar'):
   """
   Calculate the Lyman-continuum photons per second.
   """
   sp = S.FileSpectrum(specfile)  # default flux unit is flam
   if unit == "Lsolar":
      # convert the flux to erg/s/cm^2/A
      sp = sp * normfactor  # convert to erg/s/cm^2/A
   sp.convert('fnu')  # easier to work with
   freq = 3.0e8 / (sp.wave / 1.e10)  # frequency, in Hz
   n_lyc = (1. / h) * area_10pc * (sp.flux / freq)  
   # number of LyC photons as a function of wavelength
   window = (sp.wave < 912)
   NLyc = integrate.simps(n_lyc[window][::-1], x=freq[window][::-1])
   # Only integrates through wavelength < 912 A
   return np.log10(NLyc)

def add_emission_line(spec, central_wave, flux, width=20., dw=1.0):
   """
   Add a given emission line to a model spectrum, assuming Gaussian shape.
   spec is a pysynphot spectrum instance, with flux unit flam (erg/s/cm^2/A).
   Flux should therefore be given in erg/s/cm^2.
   """
   w0 = central_wave - width * 3
   w1 = central_wave + width * 3
   wave = np.arange(w0, w1+dw, dw)  # the wavelength grid for the line
   # calculate the amplitude
   # amp = flux / (np.sqrt(2.*np.pi) * width)
   line = flux * gauss.gauss(wave, central_wave, width, xstep=dw)
   line[0] = 0.; line[-1] = 0.
   new_wave = np.union1d(spec.wave, wave)  # already sorted
   spr = spec.resample(new_wave)
   sp_line = S.ArraySpectrum(wave=wave, flux=line, fluxunits='flam')
   sp_new = spr + sp_line
   return sp_new

def calc_nebular_continuum(logNLyc):
   """
   Calculate the continuum nebular emission component in erg/s/cm^2/A
   Use equations (7) and (8) from Krueger et al. 1995
   """
   c = 3.e18  # speed of light in A/s
   # free-bound to HI
   wave_HI = np.array(gamma_HI[0])
   flux_HI = (c/wave_HI**2) * (np.array(gamma_HI[1])*1.e-40/2.575e-13) * 10.**logNLyc
   neb_HI = S.ArraySpectrum(wave=wave_HI, flux=flux_HI/area_10pc, 
                            fluxunits='flam')
   # free-bound to HeI
   wave_HeI = np.array(gamma_HeI[0])
   flux_HeI = (c/wave_HeI**2) * (np.array(gamma_HeI[1])*1.e-40/2.575e-13) * 10.**logNLyc
   neb_HeI = S.ArraySpectrum(wave=wave_HeI, flux=flux_HeI*n_HeI/area_10pc,
                             fluxunits='flam')
   # Two-photon process
   wave_2q = np.array(gamma_2q[0])
   flux_2q = (c/wave_2q**2) * (np.array(gamma_2q[1])*1.e-40/2.575e-13) * 10.**logNLyc
   neb_2q = S.ArraySpectrum(wave=wave_2q, flux=flux_2q/area_10pc,
                            fluxunits='flam')
   neb_cont = neb_HI + neb_HeI + neb_2q
   return neb_cont

def add_nebular_lines(specfile, outputfile="", fluxunit='Lsolar', width=10., dw=1.0, metalKey='Z3', clobber=False, verbose=False, write_header=False, continuum=True, ebmv=0., extlaw='xgal'):
   """
   Add nebular emission lines to the model SED. The line strengths are 
   proportional to the Lyman continuum flux.
   Also can specify a value for E(B-V); this is the value for **STELLAR
   CONTINUUM**, and the value for nebular emission will be ~2x higher 
   (see Calzetti+2000). Then I attenuate both the stellar and nebular
   components before adding them up.
   """
   sp = S.FileSpectrum(specfile)
   sp = S.ArraySpectrum(sp.wave, sp.flux * (L_solar / area_10pc), 
                        fluxunits='flam')
   sp0 = S.FileSpectrum(specfile) * (L_solar / area_10pc) 
   # the reference for continuum
   # First add hydrogen recombination lines
   # Calculate LyC photon flux
   N_Lyc = 10.**(calc_Lyc(specfile, unit=fluxunit))
   # Now calculate the nebular emission component
   sp_neb = S.ArraySpectrum(sp.wave, flux=np.zeros(len(sp.wave)),
                            fluxunits='flam')
   for i in range(len(LyC_H)):
      lam = LyC_H[i][0]
      flux = LyC_H[i][1] * N_Lyc / area_10pc  # in erg/s/cm^2
      if i == 0:
         flux_Hbeta = flux
      # print lam, flux
      sp_neb = add_emission_line(sp_neb, lam, flux, width=width, dw=dw)
      # re-create an ArraySpectrum instance in order to use the resample method
      sp_neb = S.ArraySpectrum(sp_neb.wave, sp_neb.flux, fluxunits='flam')
      if verbose:
         continuum = (sp0.sample(lam-100.) + sp0.sample(lam+100.)) / 2.
         EW = flux / continuum
         print "Equiv. width at %.1f A is %.2f A." % (lam, EW)
   # Now add non-hydrogen (metal) nebular lines
   for j in range(len(nebular_ratios['lambda'])):
      lam = nebular_ratios['lambda'][j]
      flux = nebular_ratios[metalKey][j] * flux_Hbeta
      sp_neb = add_emission_line(sp_neb, lam, flux, width=width, dw=dw)
      sp_neb = S.ArraySpectrum(sp_neb.wave, sp_neb.flux, fluxunits='flam')
      if verbose:
         continuum = (sp0.sample(lam-100.) + sp0.sample(lam+100.)) / 2.
         EW = flux / continuum
         print "Equiv. width at %.1f A is %.2f A." % (lam, EW)
   if continuum:
      # Now also add nebular continuum
      neb_cont = calc_nebular_continuum(np.log10(N_Lyc))
      sp_neb = sp_neb + neb_cont
      sp_neb = S.ArraySpectrum(sp_neb.wave, sp_neb.flux, fluxunits='flam')
   sp = sp + sp_neb
   ## Now apply dust attenuation if desired (if ebmv > 0)
   if ebmv > 0:
      dust_stellar = S.Extinction(ebmv, extlaw)  # default is Calzetti law
      ebmv_neb = ebmv / 0.44  # see Calzetti+2000, Equation (3)
      dust_nebular = S.Extinction(ebmv_neb, extlaw)
      sp_stellar = sp * dust_stellar
      sp_neb = sp_neb * dust_nebular
      sp = sp_stellar + sp_neb
   if len(outputfile):
      header = ""
      output_text = ""
      if write_header:
         f = open(specfile, 'rb')
         lines = f.readlines()
         for l in lines:
            # Only keep comment lines
            if l[0] == '#':
               header += l
         f.close()
      np.savetxt(outputfile, np.array([sp.wave, sp.flux]).T, fmt='%.4f %e',
                 delimiter='   ', header=header)
      # f = open(outputfile, 'wb')
      # if write_header:
      #    for i in range(len(header)):
      #       f.write(header[i])
      # for j in range(len(sp.wave)):
      #    output_text += '%.4f   %.10f\n' % (sp.wave[j], sp.flux[j])
      #    # write flux in erg/s/cm^2/A
      #    # f.write('%.4f   %.10f\n' % (sp.wave[j], sp.flux[j]))
      # f.write(output_text)
      # f.close()
   return sp

def generate_templates(specfiles, fluxunit='Lsolar', ebmv=[]):
   """
   Generate a library of templates. Add nebular emission and add dust if 
   necessary.
   """
   # Now for each template, add dust to them
   # specfiles = glob.glob('*.sed')
   for f in specfiles:
      print f
      if not len(ebmv):
         sp = add_nebular_lines(f, f, fluxunit='Lsolar', ebmv=0)
      else:
         for i in range(len(ebmv)):
            outputfile = os.path.splitext(f)[0] + '_dust%.2f.sed' % ebmv[i]
            print outputfile
            sp = add_nebular_lines(f, outputfile, fluxunit='Lsolar', ebmv=ebmv[i])
   print "Done."

def read_age(filename):
   """
   Read the log(age) from the file name.
   """
   if 'age' not in filename:
      print "String 'age' not in file name."
      return 0
   root = os.path.splitext(filename)[0]
   log_age = root.split('age')[1]
   log_age = log_age.split('_')[0]
   return float(log_age)

def read_tau(filename):
   """
   Read the tau from the file name.
   """
   if 'tau' not in filename:
      print "String 'tau' not in file name."
      return 0
   root = os.path.splitext(filename)[0]
   tau = root.split('tau')[1]
   tau = tau.split('_')[0]
   return float(tau)

def calc_phys_prop(specfile):
   # Calculate, for each model spectrum, the physical properties to be used
   # in LePhare
   # Manual input tau (should be improved later!!)
   sp = S.FileSpectrum(specfile) * (L_solar / area_10pc)
   f = interpolate.interp1d(sp.wave, sp.flux)  # interpolated spectrum
   # c = 3.e18  angstrom/sec
   LUV = integrate.romberg(f, 2100., 2500.) / 400. * (2300.**2 / 3.e18) * area_10pc
   LR = integrate.romberg(f, 5500., 6500.) / 1000. * (6000.**2 / 3.e18) * area_10pc
   LK = integrate.romberg(f, 21000., 23000.) / 2000. * (22000.**2 / 3.e18) * area_10pc
   # LIR = -99.0
   D4000 = integrate.romberg(f, 4050, 4250) / integrate.romberg(f, 3750, 3950)
   return [LUV, LR, LK, D4000]
   # return [age, LUV, LR, LK, LIR, SMass, SFR, Z, tau, D4000]

def write_phys_prop(specfiles, physprop_file, physfile, tau_default=1.e9, Z=0.004):
   """
   To make a *.phys file for Le Phare.
   """
   # Manual input tau (should be improved later!!)
   s2 = sextractor(physprop_file)
   log_age_list = np.array(["%.4f" % (round(age, 4)) for age in s2._1])
   # print log_age_list
   ff = open(physfile, 'wb')  # the *.phys file
   i = 1  # line running count; 1st column of *.phys file
   j = 1  # model running count; 2nd column of *.phys file
   ## Le Phare seems to require a dummy line for each model with negative 
   ## values? Don't understand why...
   for f in specfiles:
      LUV, LR, LK, D4000 = calc_phys_prop(f)
      # assume the specfiles are named such that one can read the log(age) 
      # from their filenames
      log_age = f.split('age')[1][:-4]
      if 'dust' in log_age:
         log_age = log_age.split('_')[0]
      # k = np.arange(len(s2))[s2.filename==f][0]  # matching index in statFile
      k = np.arange(len(s2))[log_age_list==log_age][0]
      age = 10.**float(log_age)
      SMass = s2._7[k]
      SFR = s2._10[k]
      tau = read_tau(f)
      if tau == 0:
         tau = tau_default
      # first, write a dummy line
      ff.write("%d  %d " % (i, j))
      ff.write("%.6E  %.6E  %.6E  " % (-1, -99, -99))
      ff.write("%.6E  %.6E  %.6E  " % (-99, -99, -99))
      ff.write("%.6E  %.6E  %.6E  " % (-99, -99, tau))
      ff.write("%.6E  " % -99.0)
      ff.write('\n')
      i += 1
      # Now write the real stuff
      ff.write("%d  %d  " % (i, j))
      ff.write("%.6E  %.6E  %.6E  " % (age, np.log10(LUV), np.log10(LR)))
      ff.write("%.6E  %.6E  %.6E  " % (np.log10(LK), -99.0, SMass))
      ff.write("%.6E  %.6E  %.6E  " % (SFR, Z, tau))
      ff.write("%.6E  " % D4000)
      ff.write("\n")
      i += 1
      j += 1
   ff.close()

def merge_physfiles(physfiles, outputfile):
   """
   Merge *.phys read from different BC03 models. Need to take care of the 
   indexing to match those in the template list file in Le Phare.
   """
   mod_count = 1
   line_count = 1
   lines_all = []
   for phys in physfiles:
      with open(phys, 'rb') as f:
         lines = f.readlines()
         for i in range(len(lines)):
            l = lines[i].split()
            l = [str(line_count), str(mod_count)] + l[2:]
            lines_all += ["  ".join(l) + "\n"]
            if line_count % 2 == 0:
               mod_count += 1
            line_count += 1
   # Now write to output
   with open(outputfile, 'wb') as fout:
      for j in range(len(lines_all)):
         fout.write(lines_all[j])
   print "Done."

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

def genmod(cspname, age_min, age_max, wavemin=91.0, wavemax=1.6e6):
   """
   Automatically extract BC03 model SEDs at the age steps specified in the 
   *.4color file. Give a range of desired ages (in Gyr).
   """
   propfile = os.path.splitext(cspname)[0] + '.4color'
   props = sextractor(propfile)
   log_age = props._1
   smass = props._7
   sfr = props._10
   select = np.zeros(len(props), 'bool')
   cmd = os.getenv('bc03') + '/galaxevpl'
   f = subprocess.Popen([cmd], stdin=subprocess.PIPE)
   filenames = []
   for i in range(len(log_age)):
      age = 10.**(log_age[i] - 9.)  # in Gyr
      if (age >= age_min) & (age < age_max):
         select[i] = True
         print >> f.stdin, "%s" % cspname
         # print >> f.stdin, '%f, %f' % (wavemin, wavemax)   # enter the desired wavelength range
         print >> f.stdin, ""  # use the default value
         print >> f.stdin, "%f" % age   # should be in Gyr
         sedname = os.path.splitext(cspname)[0] + '_age%.4f.sed' % log_age[i]
         print >> f.stdin, "%s" % sedname
         filenames += [sedname]
         # retcode = f.wait()
   filenames = np.array(filenames)
   print f.poll()
   return filenames, select
   # return [filenames, log_age.take(select), smass.take(select), sfr.take(select)]
   #print retcode

def write_eazy_templist(specfiles, log_age, templist, header="", abspath=True):
   assert len(specfiles) == len(log_age)
   N = len(specfiles)
   if abspath == True:
      curdir = os.getcwd()
      specfiles = [os.path.join(curdir,x) for x in specfiles]
   count = np.arange(1, len(specfiles) + 1)
   ones = np.ones(len(specfiles))
   np.savetxt(templist,
      np.array([count,specfiles,ones,10.**(log_age-9.),ones]).T, 
      fmt=['%s','%s','%s','%s','%s'], delimiter='  ', header=header)

def write_eazy_physfile(specfiles, statfile, outputfile):
   """
   Write a *.phys file to find the best-fit physical properties for EAZY.
   statfile should be the *.4color file output from BC03.
   If I'm combining templates from different *.ised files, I'll have to merge
   the *.phys files from each of them somehow...
   """
   s2 = sextractor(statfile)
   log_age_list = np.array(["%.4f" % a for a in s2._1])
   with open(outputfile, 'wb') as f:
      f.write('# 1 FILENAME\n')
      f.write('# 2 LOG_AGE\n')
      f.write('# 3 SMASS\n')
      f.write('# 4 SFR\n')
      f.write('# 5 TAU [Gyr]\n')
      f.write('# 6 EBMV\n')
      for s in specfiles:
         root = os.path.splitext(s)[0]
         x = root.split('_')
         tau = 1.0  # default value if can't read from the file name
         ebmv = 0.0 # default value
         # First find out the age
         for y in x:
            if 'age' in y.lower():
               log_age = float(y[3:])
               break
         # Find the model index
         k = np.arange(len(s2))[log_age_list==('%.4f'%log_age)][0]
         smass = s2._7[k]
         sfr = s2._10[k]
         # Now try to read tau and ebmv from file name
         for y in x:
            if 'tau' in y:
               tau = float(y[3:])
            elif 'dust' in y:
               ebmv = float(y[4:])
         f.write('%s  %f  %e  %e  %f  %.2f\n' % (s, log_age, smass, sfr, tau, ebmv))
   print "Done."


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

