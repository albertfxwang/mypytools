#!/usr/bin/env python

# Spawned from the HSTcatalog class originally from hst_pipeline.py

import numpy as np
from pygoods import fitstable, Ftable, sextractor, angsep
from datetime import datetime
from PhotomTools import photomUtils as pu
from stats import stats_simple as ss
import matplotlib.pyplot as plt

# because Le Phare requires the input catalog contain magnitudes for each 
# filter used to create the libraries, and I don't want to create libraries 
# for each cluster, here I'll just use a union of all possible filters for 
# the SURFS UP clusters.
lephare_filters = ['f435w','f475w','f555w','f606w','f625w','f775w','f814w',\
                   'f850lp','f098m','f105w','f110w','f125w','f140w','f160w',\
                   'irac1','irac2']
# For catalogs that include IRAC fluxes, manually edit the files to add IRAC
# magnitudes/limits.
clash_bands = ['f435w','f475w', 'f555w', 'f606w','f625w','f775w','f814w','f850lp','f098m','f105w','f110w','f125w','f140w','f160w']
clash_eazy_hdr = "# id f_irac2 e_irac2 f_irac1 e_irac1 f_f160w e_f160w f_f140w e_f140w f_f125w e_f125w f_f110w e_f110w f_f105w e_f105w f_f098m e_f098m f_f850lp e_f850lp f_f814w e_f814w f_f775w e_f775w f_f625w e_f625w f_f606w e_f606w f_f555w e_f555w f_f475w e_f475w f_f435w e_f435w"
# Taken as either the central or the pivot wavelength (A)... doesn't need to 
# be too accurate because it's for plotting purposes
efflam = dict(f435w=4311., f475w=4775.7, f555w=5356.0, f606w=5888., 
              f625w=6295.5, f775w=7665.1, f814w=8115.3, f850lp=9145.2,
              f098m=9864.0, f105w=10552.0, f110w=11534.0, f125w=12486.0, 
              f140w=13923.0, f160w=15369.0, irac1=36000., irac2=45000.)


class HSTcatalog(fitstable.Ftable):
   def __init__(self, filename):
      fitstable.Ftable.__init__(self, filename)
 
   def print_sn(self, number, bands=clash_bands, fluxcol='iso'):
      # get the S/N (using scaled flux errors)
      for b in bands:
         try:
            flux = getattr(self, '%s_flux_%s' % (b.lower(), fluxcol.lower()))
            fluxerr = getattr(self, '%s_fluxerr_%s_scaled' % (b.lower(), fluxcol.lower()))
            sn = (flux / fluxerr)[self.number==number][0]
            print "%s S/N = %.2f" % (b.upper(), sn)
         except:
            print "Cannot calculate S/N in %s" % b.upper()
 
   def print_radec(self, number, print_it=False, format='degree', delimiter=':'):
      ra = self.alpha_j2000[self.number==number][0]
      dec = self.delta_j2000[self.number==number][0]
      if print_it:
         print ra, dec
      if format == 'degree':
         return ra, dec
      else:
         ra = coordinates.Angle(ra, unit='degree').hms
         ra_str = delimiter.join(['%02d'%int(ra[0]),'%02d'%int(ra[1]),'%02.2f'%ra[2]])
         dec = coordinates.Angle(dec, unit='degree').dms
         SIGN = ""
         if (dec[0] < 0) or ((dec[0]==0) and dec[1]<0):
            SIGN = "-"
         dec_str = delimiter.join([SIGN+'%02d'%abs(int(dec[0])),'%02d'%abs(int(dec[1])),'%02.2f'%abs(dec[2])])
         return [ra_str, dec_str]

   def match_radec(self, ra, dec, tol=0.5, ra_col='alpha_j2000', dec_col='delta_j2000'):
      """
      Given values of RA, DEC, find the closest object within some match
      radius (tol) in arcsec. Returns the SExtractor ID number.
      """
      racol = getattr(self, ra_col)
      deccol = getattr(self, dec_col)
      angdist = angsep.angsep(ra, dec, racol, deccol)
      if angdist.min() * 3600. > tol:
         return -1
      else:
         j = np.argsort(angdist)[0]
         return self.number[j]

   def match_radec_many(self, ra_list, dec_list, tol=0.5, ra_col='alpha_j2000', dec_col='delta_j2000'):
      """
      Return a list of SExtractor IDs for given lists of RA & DEC.
      """
      matched_id = []
      for i in range(len(ra_list)):
         matched_id.append(self.match_radec(ra_list[i], dec_list[i], tol=tol,
                           ra_col=ra_col, dec_col=dec_col))
      return matched_id

   def print_apcor(self, number, colorcol='iso', refband='f160w', print_it=True):
      # print the aperture correction from the mag. form for colors (default=ISO)
      # to the AUTO mag., using the reference band (refband)
      mag_auto = getattr(self, '%s_mag_auto' % refband)[self.number==number][0]
      mag_color = getattr(self, '%s_mag_%s' % (refband, colorcol.lower()))[self.number==number][0]
      if print_it: print mag_auto - mag_color
      return mag_auto - mag_color
 
   def get_mag(self, number, band, magform='auto', refband='f160w'):
      # also calculates aperture correction errors
      # If magform=='iso', it returns the **aperture-corrected** MAG_ISO; i.e.,
      # it will return MAG_ISO + apcor. So if band==refband, it will return the 
      # same magnitudes for magform='auto' or 'iso', but the magntiude errors
      # will be different
      band = band.lower()
      magform = magform.lower()
      mag = getattr(self, '%s_mag_%s' % (band, magform))[self.number==number][0]
      if magform == 'auto':
         apcor = 0.
         magerr_auto = 0.
      else:
         apcor = self.print_apcor(number, colorcol=magform, print_it=False,
                                  refband=refband)
         magerr_auto = getattr(self, '%s_magerr_auto' % band)[self.number==number][0]
      if magform == 'iso':
         mag_err = getattr(self, '%s_magerr_iso' % (band))[self.number==number][0]
         mag_err = np.sqrt(mag_err**2 + magerr_auto**2)
      else:
         flux = getattr(self, '%s_flux_%s' % (band,magform))[self.number==number][0]
         fluxerr = getattr(self, '%s_fluxerr_%s_scaled' % (band,magform))[self.number==number][0]
         sn = flux / fluxerr
         if sn > 0:
            mag_err = float(pu.sn2magerr(sn))
            mag_err = np.sqrt(mag_err**2 + magerr_auto**2)
         else:
            mag_err = 99.0
      try:
         mag_1sig = getattr(self, '%s_mag_%s_1sig' % (band, magform))[self.number==number][0]
      except:
         mag_1sig = mag_err
      if magform.lower() != 'auto':
         mag = mag + apcor
      if (mag >= 90.0) | (mag_err>=1.0857):
         return 99.0, mag_1sig + apcor
      else:
         return mag, mag_err

   def get_all_mags(self, number, bands, colorcol='iso'):
      # Assume that bands[0] is the reference band
      mags = np.zeros(len(bands))
      magerrs = np.zeros(len(bands))
      for i in range(len(bands)):
         if i == 0:
            mag, magerr = self.get_mag(number, bands[i], magform='auto', 
                                       refband=bands[0])
         else:
            mag, magerr = self.get_mag(number, bands[i], magform=colorcol,
                                       refband=bands[0])
         mags[i] = mag
         magerrs[i] = magerr
      return mags, magerrs
 
   def print_mag(self, number, band, magform='auto', print_it=True, refband='f160w'):
      # Does NOT include errors in aperture correction
      mag, mag_err = self.get_mag(number, band, magform=magform, 
                                  refband=refband)
      outstr = "%.2f +/- %.3f" % (mag, mag_err)
      if print_it:
        print outstr
      return outstr 

   def plot_all_mags(self, objid, bands, colorcol='iso', ax=None, xoffset=0., **ebar_kwargs):
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      mags, magerrs = self.get_all_mags(objid, bands, colorcol=colorcol)
      wave = np.array([efflam[b] for b in bands]) + xoffset
      ax.errorbar(wave, mags, yerr=magerrs, **ebar_kwargs)
      plt.gca().invert_yaxis()
      return ax

   def plot_all_mags_objects(self, objids, bands, objnames=None, colorcol='iso', **ebar_kwargs):
      ebar_kwargs_default = dict(elinewidth=2, capsize=8)
      for k in ebar_kwargs_default.keys():
         if k not in ebar_kwargs.keys():
            ebar_kwargs[k] = ebar_kwargs_default[k]
      fig = plt.figure()
      ax = fig.add_subplot(111)
      for i in range(len(objids)):
         x = objids[i]
         xoffset = 20 * i
         if objnames == None:
            objname = '%d' % x
         else:
            assert len(objnames) == len(objids), "Length of objnames is different from the length of objids..."
            objname = objnames[i]
         ax = self.plot_all_mags(x, bands, colorcol=colorcol, ax=ax, 
                                 xoffset=xoffset, label=objname, **ebar_kwargs)
      ax.legend(loc=0)
      ax.set_xlabel('Wavelength (A)')
      ax.set_ylabel('magnitude [%s]' % colorcol.upper())
      return ax

   def combine_mag(self, numbers, band, magzero, magform='auto', print_it=False):
      # Combine fluxes from more than 1 segmentation
      nobj = len(numbers)
      fluxes = np.zeros(nobj)
      fluxerrs = np.zeros(nobj)
      for i in range(nobj):
         # Read the fluxes
         num = numbers[i]
         fluxes[i] = getattr(self,'%s_flux_%s'%(band,magform))[self.number==num][0]
         fluxerrs[i] = getattr(self,'%s_fluxerr_%s'%(band,magform))[self.number==num][0]
         apcor = self.print_apcor(num, colorcol=magform, print_it=False)
         fluxcor = 10.**(-0.4 * apcor)
         if (fluxes[i] / fluxerrs[i]) > 1:
            fluxes[i] = fluxes[i] * fluxcor
         else:
            fluxerrs[i] = fluxerrs[i] * fluxcor
      # Now combine the flux and flux error
      fluxTot = fluxes.sum()
      fluxerrTot = np.sqrt((fluxerrs**2).sum())
      if (fluxTot / fluxerrTot) > 1.0:
         magTot = magzero - 2.5 * np.log10(fluxTot)
         magerrTot = float(pu.sn2magerr(fluxTot / fluxerrTot))
      else:
         magTot = 99.0
         magerrTot = magzero - 2.5 * np.log10(fluxerrTot)
      return magTot, magerrTot

   def calc_color(self, objid, band1, band2, colorcol='iso', print_it=False):
      mag1, magerr1 = self.get_mag(objid, band1, magform=colorcol)
      mag2, magerr2 = self.get_mag(objid, band2, magform=colorcol)
      if (mag1 < 90) & (mag2 < 90):
         color = mag1 - mag2
         color_err = np.sqrt(magerr1**2 + magerr2**2)
      else:
         color = 0.
         color_err = 10.
      return color, color_err

   def calc_all_colors(self, objid, bands, colorcol='iso', print_it=True):
      colors = []
      color_errs = []
      for i in range(len(bands) - 1):
         color, color_err = self.calc_color(objid, bands[i], bands[i+1],
                                            colorcol=colorcol)
         colors.append(color)
         color_errs.append(color_err)
         if print_it:
            print "%s-%s = %.3f +/- %.3f" % (bands[i],bands[i+1],color,color_err)
      return np.array(colors), np.array(color_errs)

   def plot_all_colors(self, objid, bands, ax=None, colorcol='iso', xoffset=0., **ebar_kwargs):
      """
      Plot all colors for one object.
      """
      wave = np.array([efflam[b] for b in bands])
      # Use the mean wavelength between adjacent filters for the color
      # colorwave = np.array([np.mean([wave[i],wave[i+1]]) for i in range(len(wave)-1)]) + xoffset
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      colors, color_errs = self.calc_all_colors(objid, bands, colorcol=colorcol)
      xdata = np.arange(len(wave)-1) + 0.5 + xoffset
      ax.errorbar(xdata[color_errs<10], colors[color_errs<10], 
                  yerr=color_errs[color_errs<10], **ebar_kwargs)
      # set x-axis labels
      ax.set_xticks(range(len(wave)))
      ax.set_xticklabels([b.upper() for b in bands])
      ax.set_xlabel('Filters')
      ax.set_xlim(xmin=-0.5, xmax=len(wave)-0.5)
      ax.set_ylabel('Color [%s]' % colorcol.upper())
      return ax

   def plot_all_colors_objects(self, objids, bands, objnames=None, colorcol='iso', **ebar_kwargs):
      ebar_kwargs_default = dict(elinewidth=2, capsize=8, linestyle='none',
                                 fmt='x', ms=10, mew=2)
      for k in ebar_kwargs_default.keys():
         if k not in ebar_kwargs.keys():
            ebar_kwargs[k] = ebar_kwargs_default[k]
      fig = plt.figure()
      ax = fig.add_subplot(111)
      for i in range(len(objids)):
         x = objids[i]
         xoffset = 0.05 * i
         if objnames == None:
            objname = '%d' % x
         else:
            assert len(objnames) == len(objids), "Length of objnames is different from the length of objids..."
            objname = objnames[i]
         ax = self.plot_all_colors(x, bands, colorcol=colorcol, ax=ax, 
                                 xoffset=xoffset, label=objname, **ebar_kwargs)
      ax.legend(loc=0)
      ax.set_ylabel('Color [%s]' % colorcol.upper())
      return ax


   def test_same_color(self, objids, bands, colorcol='iso'):
      """
      Test if all the colors from objids are consistent within 
      errors.
      """
      print "Comparing colors among %s:" % str(objids)
      for i in range(len(bands) - 1):
         b1 = bands[i]
         b2 = bands[i+1]
         colors = [self.calc_color(x, b1, b2, colorcol=colorcol)[0] for x in objids]
         color_errs = [self.calc_color(x, b1, b2, colorcol=colorcol)[1] for x in objids]
         # color1, color_err1 = self.calc_color(objid1, b1, b2)
         # color2, color_err2 = self.calc_color(objid2, b1, b2)
         # same = ss.consistent_sigma(color1, color_err1, color2, color_err2)
         same = ss.consistent_sigma(colors, color_errs)
         print "%s-%s are consistent? %s." % (b1, b2, str(same))

 
   def print_mag_objects(self, numbers, band, magform='auto'):
      # DOES include errors in aperture correction if magform != 'auto'
      if not hasattr(self, '%s_mag_auto' % band):
         print "Magnitudes for %s does not exist." % band.upper()
         return 0
      if magform == 'auto':
         mags = map(lambda x: self.get_mag(x, band)[0], numbers)
         magerrs = map(lambda x: self.get_mag(x, band)[1], numbers)
      else:
         mags = map(lambda x: self.get_mag(x, band, magform=magform, )[0], numbers)
         magerrs = map(lambda x: self.get_mag(x, band, magform=magform)[1], numbers)
      mags_str = ','.join(map(lambda y: '%.4f'%y, mags))
      magerrs_str = ','.join(map(lambda y: '%.4f'%y, magerrs))
      print "%s magnitudes:" % band.upper()
      print mags_str
      print "%s magnitude errors:" % band.upper()
      print magerrs_str
 
   def print_radec_objects(self, numbers):
      ra = map(lambda x: self.print_radec(x, print_it=False)[0], numbers)
      dec = map(lambda x: self.print_radec(x, print_it=False)[1], numbers)
      ra_str = ','.join(map(lambda y: str(y), ra))
      dec_str = ','.join(map(lambda y: str(y), dec))
      print "RA:"
      print ra_str
      print "DEC:"
      print dec_str
 
   def print_all(self, objid, bands, colorcol='iso'):
      # print the following information about the source with objid
      # RA, DEC
      # MAG_AUTO from the reference band (refband; the FIRST element in bands)
      # MAG_ISO from other bands, normalized to match MAG_AUTO in refband 
      print "Information about object %d:" % objid
      print "RA, DEC = %.7f, %.7f" % self.print_radec(objid, print_it=False) 
      mags = []
      magerrs = []
      apcor = self.print_apcor(objid, colorcol=colorcol, refband=bands[0], print_it=False)
      print "%s magnitude: %s" % (bands[0], self.print_mag(objid, bands[0], magform='auto', print_it=False))
      for b in bands[1:]:
         try:
            print "%s magnitude: %s" % (b, self.print_mag(objid, b, magform=colorcol, 
                                      print_it=False))
         except:
            print "Could not get magnitude from %s..." % b.upper()
 
   def write_reg_column(self, filename, column, radius=1.0, color='green'):
      # write a region file (use circular regions; radius is in arcsec)
      f = open(filename, 'wb')
      col = getattr(self, column)
      f.write('global color=%s width=2\n' % color)
      for i in range(len(col)):
         if col[i]:
            f.write('fk5; circle(%f, %f, %.2f") # text={%d}\n' % (self.alpha_j2000[i], self.delta_j2000[i], radius, self.number[i]))
      f.close()
 
   def write_reg_obj(self, filename, objid, radius=2.0, color='green'):
      f = open(filename, 'wb')
      f.write('global color=%s\n' % color)
      for x in objid:
         ra, dec = self.print_radec(x, print_it=False)
         f.write('fk5; circle(%f, %f, %.2f") # text={%d}\n' % (ra, dec, radius, x))
      f.close()
 
   def print_eazy_flux(self, objid, eazyid=[], header=clash_eazy_hdr, detect_band='f160w', SNLim=1.0, colorcol='iso'):
      # Provide a header line for the flux catalog to be fed into EAZY, print 
      # the fluxes and errors in uJy for the objects specified in objid
      # First, parse the header line; assume it starts with #
      # First column is always the ID; then the subsequent columns are in the 
      # format of f_[filter_name] and e_[filter_name]
      # make sure there is a zphot.translate to convert the filter names into
      # something EAZY can understand.
      h = header.split()
      obj_str = ""
      if len(eazyid):
         assert len(objid)==len(eazyid), "Please provide the same number for EAZY ID as the SExtractor ID."
      for i in range(len(objid)):
         if len(eazyid):
            obj_str += "%s  " % eazyid[i]
         else:
            obj_str += "%d  " % objid[i]
         for j in range(len(h))[2:]:
            # apcor = self.print_apcor(objid[i], print_it=False)
            if h[j].startswith('f_'):
               b = h[j][2:].lower()
               if b == detect_band:
                  mag, mag_err = self.get_mag(objid[i], b, magform='auto')
               else:
                  try:
                     mag, mag_err = self.get_mag(objid[i], b, magform=colorcol)
                  except:
                     # raise KeyError, "Filter %s does not exist in the catalog." % b
                     mag, mag_err = 99.0, -99.0
               # Now format the magnitudes to comply with EAZY convention
               if mag > 90.:  # not detected
                  if mag_err > 0:  # there is flux limit; mag_err gives the 1-sigma limit
                     flux = 0.
                     flux_err = pu.ABmag2uJy(mag_err) * SNLim  # report N-sigma flux limit
                  else:
                     flux = -99.
                     flux_err = -99.
               else:
                  # Test if S/N is higher than SNLim
                  SN = pu.magerr2sn(mag_err)
                  if SN >= SNLim:
                     flux = pu.ABmag2uJy(mag)
                     sn = pu.magerr2sn(mag_err)
                     flux_err = flux / sn
                  else:
                     # convert to the N-sigma upper limit
                     flux = pu.ABmag2uJy(mag)  # calculate the reported flux first
                     flux_1sig = flux / SN  # calculate 1-sigma flux limit
                     flux = 0.  # set flux to zero to enforce flux upper limit
                     flux_err = flux_1sig * SNLim  # calculate N-sigma upper limit
               obj_str += "%.6e  %.6e  " % (flux, flux_err)
            else:
               pass
         obj_str += "\n"
      return obj_str
 
   def irac_mag_to_uJy(self, mag, magerr, print_it=True):
     if mag < 90:
       flux = pu.ABmag2uJy(mag)
       sn = pu.magerr2sn(magerr)
       fluxerr = flux / sn
     else:
       if magerr < 0:
         flux = -99.0
         fluxerr = -99.0
       else:
         flux = 0.
         fluxerr = pu.ABmag2uJy(magerr)
     if print_it:
        print "%.6e  %.6e" % (flux, fluxerr)
     return flux, fluxerr

   def irac_mag_to_uJy_2bands(self, mag2, mag2err, mag1, mag1err):
      flux2, flux2err = self.irac_mag_to_uJy(mag2, mag2err, print_it=False)
      flux1, flux1err = self.irac_mag_to_uJy(mag1, mag1err, print_it=False)
      return "%.6e  %.6e  %.6e  %.6e" % (flux2, flux2err, flux1, flux1err)
 
   def write_LePhare_flux(self, objid, filterNames=lephare_filters, refBand='f160w', output="", SNLim=1.0):
      """
      Write a photometric catalog in the format accepted by Le Phare photo-z code.
      Works for CAT_TYPE = SHORT. Assume CAT_MAG = AB and CAT_FMT = MEME.
      """
      hdrString = "# Current time: " + str(datetime.now()) + "\n"
      hdrString += "# ID "
      for f in filterNames:
        hdrString += "%s_mag %s_magerr " % (f.lower(), f.lower())
      hdrString += "\n"
      objString = ""
      for i in range(len(objid)):
         objString += "%d  " % objid[i]
         for j in range(len(filterNames)):
            fname = filterNames[j].lower()
            try:
               # apcor = self.print_apcor(objid[i], refband=refBand, print_it=False)
               if fname == refBand:  # in which we calculate aperture correction
                  mag, magErr = self.get_mag(objid[i], fname, magform='auto')
               else:
                  try:
                     mag, magErr = self.get_mag(objid[i], fname, magform='iso')
                  except:
                     print "Filter %s does not exist in the catalog." % fname.upper()
                     mag, magErr = [-99., -99.]
               # Now also determine if the magnitudes and errors make sense
               # first, check if there is valid flux measurement
               if magErr <= 0.:
                  # no information for this band -- set both mag and error to -99.0
                  mag, magErr = [-99., -99.]
               elif mag > 90:
                  # not detected (mag = 99); set N-sigma magnitude limit, where
                  # N = SNLim
                  mag = magErr - 2.5 * np.log10(SNLim)
                  magErr = -1.0
               else:
                  #test if S/N>=SNLim
                  SN = pu.magerr2sn(magErr)
                  print objid[i], fname, SN, mag, magErr
                  if SN < SNLim:
                     # flux = pu.ABmag2uJy(mag)
                     # flux_1sig = flux / SN
                     # flux_nsig = flux_1sig * SNLim
                     # mag = pu.uJy2ABmag(flux_nsig)
                     mag = pu.calcNsigMag(mag, magErr, N=SNLim)
                     magErr = -1.0
                  print mag, magErr
               # Now append to objString
               # print "%s, %.3f, %.3f" % (fname, mag, magErr)
               objString += "%.4f %.4f " % (mag, magErr)
            except:
               # this object is not observed in this filter.
               objString += "-99.0  -99.0  "
         objString += "\n"
      if len(output):
         f = open(output, 'wb')
         f.write(hdrString)
         f.write(objString)
         f.close()
      else:
         print hdrString
         print objString 