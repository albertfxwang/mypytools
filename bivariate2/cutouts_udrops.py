#!/usr/bin/env python

import numpy as np
from ImageTools import cutouts
from pygoods import sextractor, Ftable

drz_u = '/Users/khuang/CANDELS/goodss/mosaics/vimos_u/tfit_015_025_sqr_6_bg1.fits'
drz_f435w = '/Users/khuang/CANDELS/goodss/mosaics/goods_s_acs_v3/gs_presm4_all_acs_f435w_60mas_v3.0_drz.fits'
drz_f606w = '/Users/khuang/CANDELS/goodss/mosaics/goods_s_acs_v3/gs_presm4_all_acs_f606w_60mas_v3.0_drz.fits'
drz_f098m = '/Users/khuang/CANDELS/goodss/mosaics/all_combined_v0.5/gs_all_candels_ers_f098m_060mas_v0.5_drz.fits'
drz_f105w = '/Users/khuang/CANDELS/goodss/mosaics/all_combined_v0.5/gs_all_candels_ers_udf_f105w_060mas_v0.5_drz.fits'
drz_f160w = '/Users/khuang/CANDELS/goodss/mosaics/all_combined_v0.5/gs_all_candels_ers_udf_f160w_v0.5_drz.fits'
seg_f160w = '/Users/khuang/CANDELS/goodss/mosaics/all_combined_v0.5/gs_all_sx_h_120604_hphotom_comb_seg_psfmatch2h.fits'
images = [drz_f435w, drz_f606w, drz_f098m, drz_f105w, drz_f160w, seg_f160w]
filters = ['f435w', 'f606w', 'f098m', 'f105w', 'f160w', 'seg_f160w']
catalog_udrops = '/Users/khuang/Dropbox/Research/bivariate/udrops_sample/gds_udrops_all_140313.fits'

class UdropsCutouts(cutouts.Cutouts):
   def __init__(self, images=images, filters=filters, catalog=catalog_udrops, format='fits', objid='objid'):
      super(UdropsCutouts, self).__init__(images, filters)
      self.use_catalog(catalog, format=format, objid=objid)

   def use_catalog(self, catalog, objid='id', ra='ra', dec='dec', format='fits'):
      if format.lower() == 'fits':
         self.c = Ftable(catalog)
         self.Nc = len(self.c.d)
      else:
         self.c = sextractor(catalog)
         self.Nc = len(self.c)
      self.objid = getattr(self.c, objid)
      self.ra = getattr(self.c, ra)
      self.dec = getattr(self.c, dec)

   def cut_objid(self, objid, width):
      # Make cutouts using object ID.
      assert objid in self.objid, "Object ID %d not found." % objid
      ra = self.ra[self.objid==objid][0]
      dec = self.dec[self.objid==objid][0]
      name = 'obj%d' % objid
      self.cut_radec_all(ra, dec, self.filters, width, name, norm=False)
