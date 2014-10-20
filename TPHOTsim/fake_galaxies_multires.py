import numpy as np
from sexsim import fake_galaxies as fg
import pywcs
import pyfits

"""
Generate fake galaxies in multi-resolution images, in preparation for 
TPHOT/TFIT simulations.
For the moment, we only insert fake galaxies into the high-res images, so the 
part that inserts fake galaxies into low-res images will be incomplete.

We assume that the sources are UNRESOLVED in the low-res images.

Created: 2014/07/15
"""

class FakeGalaxiesMultiRes(fg.FakeGalaxies):
   """
   Generate random positions based on RA, DEC instead of just image 
   coordinates. This way one can put in galaxies at consistent locations 
   between multi-resolution images.
   """
   def __init__(self, *args, **kwargs):
      fg.FakeGalaxies.__init__(self, *args, **kwargs)

   def get_radec(self, edgebuffer=60, flagmax=1):
      print "Randomly picks RA, DEC for each fake source from the high-res image..."
      # hdr = pyfits.getheader(self.realimages[self.bands[0]])
      # self.wcs = pywcs.WCS(hdr)
      hdr_hires = pyfits.getheader(self.realimages[self.bands[0]])
      hdr_lores = pyfits.getheader(self.realimages[self.bands[1]])
      self.wcs_hires = pywcs.WCS(hdr_hires)
      self.wcs_lores = pywcs.WCS(hdr_lores)
      flag_img = pyfits.getdata(self.flagimages[self.bands[0]])
      x0 = edgebuffer/2.; y0 = edgebuffer/2.
      x1 = self.xmax - edgebuffer/2.; y1 = self.ymax - edgebuffer/2.
      ra0, dec0 = self.wcs_hires.wcs_pix2sky([[x0, y0]], 1)[0]
      ra1, dec1 = self.wcs_hires.wcs_pix2sky([[x1, y1]], 1)[0]
      ra_arr = np.zeros(self.ngal)
      dec_arr = np.zeros(self.ngal)
      for i in range(self.ngal):
         offimage = 1
         while offimage:
            ra = np.random.uniform(ra0, ra1)
            dec = np.random.uniform(dec0, dec1)
            x, y = self.wcs_hires.wcs_sky2pix([[ra, dec]], 1)[0]
            if fg.test_flag(x, y, flag_img) < flagmax:
               offimage = 0
               ra_arr[i] = ra
               dec_arr[i] = dec
      self.ra = ra_arr
      self.dec = dec_arr
      return ra_arr, dec_arr

   def get_radec_around(self, RA, DEC, radius=10.0, flagmax=1, inner_radius=1.0):
      """
      Get random RA, DEC around a list of positions, one random position for 
      each object. The maximum separation from each source is specified by 
      the keyword argument radius (in arcsec).

      RA, DEC should be in DEGREES.
      """
      # generate a random position by randomly drawing a radius and a position
      # angle
      flag_img = pyfits.getdata(self.flagimages[self.bands[0]])
      radius_deg = radius / 3600.
      ra_arr = np.zeros(self.ngal)
      dec_arr = np.zeros(self.ngal)
      for i in range(self.ngal):
         offtarget = 1
         while offtarget:
            R = np.random.uniform(0., 1.) * radius_deg
            R = np.maximum(R, inner_radius / 3600.)
            THETA = np.random.uniform(0., 360.)
            ra = RA[i] - R * np.cos(THETA)
            dec = DEC[i] + R * np.sin(THETA)
            x, y = self.wcs_hires.wcs_sky2pix([[ra, dec]], 1)[0]
            if fg.test_flag(x, y, flag_img) < flagmax:
               offtarget = 0
               ra_arr[i] = ra
               dec_arr[i] = dec
      self.ra = ra_arr
      self.dec = dec_arr
      return ra_arr, dec_arr

   def get_xy(self, RA=[], DEC=[], edgebuffer=60, flagmax=1, mode='hires'):
      """
      Returns the image coordinates of the fake sources. If values for RA, DEC
      are given, use those to calculate the image coordinates. Otherwise, 
      randomly choose RA, DEC for each source and return their image 
      coordinates.

      Keyword argument "mode" specifies which image to get image coordinate 
      from (hires or lores).
      """
      assert len(RA) == len(DEC)
      if len(RA):
         radec = np.array([RA, DEC]).T
      else:
         ra_arr, dec_arr = self.get_radec(edgebuffer=edgebuffer, 
                                          flagmax=flagmax)
         radec = np.array([ra_arr, dec_arr]).T
      if mode == 'lores':
         xy = self.wcs_lores.wcs_sky2pix(radec, 1)
      else:
         xy = self.wcs_hires.wcs_sky2pix(radec, 1)
      xarr = xy[:,0]
      yarr = xy[:,1]
      
      return xarr, yarr

   def makegals_multiband(self, flagimage, igalfile=""):
      # Make artdata object list for the high-res image
      fg.FakeGalaxies.makegals_multiband(self, flagimage, igalfile=igalfile,
                                         bands=[self.bands[0]])

   def makegals_hires(self, flagimage, igalfile=""):
      self.makegals_multiband(flagimage, igalfile=igalfile)

   def makegals_lores(self):
      """Makes the galaxy list files for the low-res image. """
      
      print "in makegals_lores"      
      # Write the galaxy parameters out to files
      self.artfiles[self.bands[1]] = "glart_%s.list" % (self.bands[1])  
      # input file for iraf.mkobject
      f = open(self.artfiles[self.lores_bands[1]], "wb")
      for i in range(self.ngal):
         # Write lines for the mkobjects galaxy list
         if self.gtype[i] == devauc:
            f.write("%10.2f %10.2f %8.3f %12s %8.3f %6.2f %6.2f no " \
               % (self.x_lores[i], self.y_lores[i], \
                  self.mag[self.bands[1]][i]+devauc_correction, \
                  gtype_str(self.gtype[i]), self.re[i], self.axis_ratio[i], \
                  self.position_angle[i]))
         else:
            f.write("%10.2f %10.2f %8.3f %12s %8.3f %6.2f %6.2f no " 
               % (self.x_lores[i], self.y_lores[i],\
                  self.mag[self.bands[1]][i], gtype_str(self.gtype[i]), \
                  self.re[i], self.axis_ratio[i], self.position_angle[i]))
         f.write("\n")
      
      f.close()
      # Write out galaxies to igalfile if desired
      if len(igalfile) > 0:
         igfile = open(igalfile, 'a')  # the *.allgal file in the detection-band directory
         for i in range(self.ngal):
            outstring = "%10.2f %10.2f " % (self.x_lores[i], self.y_lores[i])
            outstring = outstring + "%d %8.3f %6.2f %6.2f" % (self.gtype[i], \
                        self.re[i], self.axis_ratio[i], self.position_angle[i])
            outstring = outstring + "%8.3f " % (self.mag[self.bands[1]][i])
            if len(self.othercols):
               keys = self.othercols.keys()
               for k in keys:
                  outstring = '%s ' % self.othercols[k][i]
            igfile.write("%s\n" % outstring)
         igfile.close()
      print "finish makegals_lores"

