#!/usr/bin/env python

# Spawned from the HSTcatalog class originally from hst_pipeline.py

import numpy as np
from pygoods import fitstable, Ftable, sextractor
from datetime import datetime
from PhotomTools import photomUtils as pu

# because Le Phare requires the input catalog contain magnitudes for each 
# filter used to create the libraries, and I don't want to create libraries 
# for each cluster, here I'll just use a union of all possible filters for 
# the SURFS UP clusters.
lephare_filters = ['f435w','f475w','f555w','f606w','f625w','f775w','f814w',\
                   'f850lp','f098m','f105w','f110w','f125w','f140w','f160w',\
                   'irac1','irac2']
# For catalogs that include IRAC fluxes, manually edit the files to add IRAC
# magnitudes/limits.
clash_bands = ['f435w','f475w', 'f555w', 'f606w','f625w','f775w','f814w','f850lp','f105w','f110w','f125w','f140w','f160w']


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
 
   def print_apcor(self, number, colorcol='iso', refband='f160w', print_it=True):
      # print the aperture correction from the mag. form for colors (default=ISO)
      # to the AUTO mag., using the reference band (refband)
      mag_auto = getattr(self, '%s_mag_auto' % refband)[self.number==number][0]
      mag_color = getattr(self, '%s_mag_%s' % (refband, colorcol.lower()))[self.number==number][0]
      if print_it: print mag_auto - mag_color
      return mag_auto - mag_color
 
   def get_mag(self, number, band, magform='auto'):
      # also calculates aperture correction errors
      band = band.lower()
      magform = magform.lower()
      mag = getattr(self, '%s_mag_%s' % (band, magform))[self.number==number][0]
      if magform == 'auto':
         apcor = 0.
         magerr_auto = 0.
      else:
         apcor = self.print_apcor(number, colorcol=magform, print_it=False)
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
 
   def print_mag(self, number, band, magform='auto', print_it=True):
      # Does NOT include errors in aperture correction
      mag, mag_err = self.get_mag(number, band, magform=magform)
      outstr = "%.2f +/- %.3f" % (mag, mag_err)
      if print_it:
        print outstr
      return outstr 

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

   def calc_color(self, objid, band1, band2, print_it=False):
      mag1, magerr1 = self.get_mag(objid, band1, magform='iso')
      mag2, magerr2 = self.get_mag(objid, band2, magform='iso')
      color = mag1 - mag2
      color_err = np.sqrt(magerr1**2 + magerr2**2)
      return color, color_err
 
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
      f.write('global color=%s\n' % color)
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
 
   def print_eazy_flux(self, header, objid, detect_band='f160w', SNLim=3.0):
      # Provide a header line for the flux catalog to be fed into EAZY, print 
      # the fluxes and errors in uJy for the objects specified in objid
      # First, parse the header line; assume it starts with #
      # First column is always the ID; then the subsequent columns are in the 
      # format of f_[filter_name] and e_[filter_name]
      # make sure there is a zphot.translate to convert the filter names into
      # something EAZY can understand.
      h = header.split()
      obj_str = ""
      for i in range(len(objid)):
         obj_str += "%d  " % objid[i]
         for j in range(len(h))[2:]:
            # apcor = self.print_apcor(objid[i], print_it=False)
            if h[j].startswith('f_'):
               b = h[j][2:].lower()
               if b == detect_band:
                  mag, mag_err = self.get_mag(objid[i], b, magform='auto')
               else:
                  try:
                     mag, mag_err = self.get_mag(objid[i], b, magform='iso')
                  except:
                     raise KeyError, "Filter %s does not exist in the catalog." % b
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
 
   def write_LePhare_flux(self, objid, filterNames=lephare_filters, refBand='f160w', output="", SNLim=3.0):
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
                  print objid[i], fname, SN
                  if SN < SNLim:
                     # flux = pu.ABmag2uJy(mag)
                     # flux_1sig = flux / SN
                     # flux_nsig = flux_1sig * SNLim
                     # mag = pu.uJy2ABmag(flux_nsig)
                     mag = pu.calcNsigMag(mag, magErr)
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