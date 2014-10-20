#!/usr/bin/env python

from numpy import *
from pygoods import *
import os, sys, glob


#def colorsel_g04(bmags,vmags,zmags):
#   """Giavalisco et al. 2004 criteria"""
#   bmv = bmags - vmags
#   vmz = vmags - zmags
#   c1 = (bmv >= (1.2 + 1.4 * vmz))
#   c2 = (bmv >= 1.2)
#   c3 = (vmz <= 1.2)
#   return (c1*c2*c3)
def colorsel_g04(bflux, bfluxerr, bzpt, vflux, vfluxerr, vzpt, zflux, zfluxerr, zzpt):
	"""Giavalisco et al. 2004 criteria"""
	c0 = (zflux / zfluxerr >= 5.0)   # z-band S/N criteria
	bmags = flux2mag_ul(bflux,bfluxerr,bzpt)
	vmags = flux2mag_ul(vflux,vfluxerr,vzpt)
	zmags = flux2mag(zflux,zzpt)
	bmv = bmags - vmags
	vmz = vmags - zmags
	c1 = (bmv >= (1.2 + 1.4 * vmz))
	c2 = (bmv >= 1.2)
	c3 = (vmz <= 1.2)
	return (c0 & c1 & c2 & c3)
	

#def colorsel_l09(bmags,vmags,zmags):
#   """S.-K. Lee et al. 2009 criteria"""
#   bmv = bmags - vmags
#   vmz = vmags - zmags
#   c1 = (bmv > (1.1 + vmz))
#   c2 = (bmv > 1.1)
#   c3 = (vmz < 1.6)
#   color_crit = (c1 * c2 * c3)
def colorsel_l09(bflux, bfluxerr, bzpt, vflux, vfluxerr, vzpt, zflux, zfluxerr, zzpt):
	"""S.-K. Lee et al. 2009 criteria"""
	c0 = (zflux / zfluxerr >= 5.0)   # z-band S/N criteria
	bmags = flux2mag_ul(bflux,bfluxerr,bzpt)
	vmags = flux2mag_ul(vflux,vfluxerr,vzpt)
	zmags = flux2mag(zflux,zzpt)
	bmv = bmags - vmags
	vmz = vmags - zmags
	c1 = (bmv > (1.1 + vmz))
	c2 = (bmv > 1.1)
	c3 = (vmz < 1.6)
	return (c0 & c1 & c2 & c3)


def colorsel_g04_exp(bflux,bfluxerr,bzpt,vflux,vfluxerr,vzpt,zflux,zfluxerr,zzpt):
   """Giavalisco et al. 2004 criteria, expanded by 0.2 mag in all directions."""
   c0 = (zflux / zfluxerr >= 5.0)   # z-band S/N criteria
   bmags = flux2mag_ul(bflux,bfluxerr,bzpt)
   vmags = flux2mag_ul(vflux,vfluxerr,vzpt)
   zmags = flux2mag(zflux,zzpt)
   bmv = bmags - vmags
   vmz = vmags - zmags
   c1 = (bmv >= (1.0 + 1.4 * vmz))
   c2 = (bmv >= 1.4)
   c3 = (vmz <= 1.0)
   return (c0 & c1 & c2 & c3)


def flux2mag_ul(flux,fluxerr,zpt):
	# convert detections with S/N < 1 sigma to upper limits
   aflux = where((flux<=0)|(flux<=fluxerr), fluxerr, flux)  # if flux<0 or S/N<1, use fluxerr
   mag = zpt - 2.5 * log10(aflux)
   return mag
   
def flux2mag(flux, zpt):
	# convert flux straight to mag, without worrying about upper limits
	# beware of negative values
	mag = zpt - 2.5 * log10(flux)
	mag = where(flux<0, 99.0, mag)
	return mag


## Below are wrapper scripts that take specific catalogs for selection... these are only
## kept below for reference. They are not maintained after 6/18/12.


def bdropsel(cb, cv, cz, cs_thresh=0.95, zmagcut=26.5, imaflags_thresh=4):
   return dropsel_magiso(cb, cv, cz, cs_thresh=cs_thresh, zmagcut=zmagcut,
      imaflags_thresh=imaflags_thresh)


def dropsel_magiso(bcat,vcat,zcat,cs_thresh=0.95,zmagcut=26.5,
   imaflags_thresh=4):
   """Select objects with non-stellar morphology, good image flags,
      z-band magnitude cuts if necessary, weed out negative fluxes,
      then do color selection
      cs_thresh: class_star threshold (only keep points below this value)
      # output a boolean array of the same length as input catalogs
      # with True values meaning it's selected as a B-dropout
   """
   # class_star threshold down to zmag = 26.0 only
   cs_crit = (zcat.class_star<=cs_thresh) | (zcat.mag_auto > 26.0)
   # z-band magnitude cut
   zmag_crit = (zcat.mag_auto<=zmagcut)
   # image flags
   flag_crit = (vcat.imaflags_iso<=imaflags_thresh)&(zcat.imaflags_iso<=imaflags_thresh)
   # S/N cut - clarification on how to calculate S/N?
   zston_crit = ((zcat.flux_auto/zcat.fluxerr_auto)>=5.)
   # color selection
   bmags = bcat.mag_iso.copy()
   vmags = vcat.mag_iso.copy()
   zmags = zcat.mag_iso.copy()
   # Use 1-sigma as upper limits in B-band if S/N <= 1.0
   for i in range(len(bcat)):
      if (bcat.flux_iso[i]/bcat.fluxerr_iso[i]) <= 1.0:
         bmags[i] = flux2mag(bcat.fluxerr_iso[i],25.65288)
      if (vcat.flux_iso[i]/vcat.fluxerr_iso[i]) <= 1.0:
         vmags[i] = flux2mag(vcat.fluxerr_iso[i],26.49341)

   color_crit = colorsel_g04(bmags,vmags,zmags)

   crit = cs_crit & zmag_crit & color_crit & zston_crit & flag_crit
   return crit


def make_bdropcat(bcat,vcat,icat,zcat,outname,photzcat,cullcat,cs_thresh=0.95,zmagcut=26.5,
   imaflags_thresh=4,magform='magiso',tol_as=0.5):
   # makes B-dropout catalog in either North or South
   # magform -- format of SE magnitude used in color selection (mag_auto or mag_iso)
   #            it's a very important choice, and still not clear for now which one should
   #            be used
   # tol_as -- matching radius in arcsec when matching with spec-z catalog in GOODS-N
   bcat = sextractor(bcat)
   vcat = sextractor(vcat)
   icat = sextractor(icat)
   zcat = sextractor(zcat)
   # do color selection
   if magform == 'magauto':
      crit = dropsel_magauto(bcat,vcat,zcat,cs_thresh=cs_thresh,zmagcut=zmagcut,
         imaflags_thresh=imaflags_thresh)
   elif magform == 'magiso':
      crit = dropsel_magiso(bcat,vcat,zcat,cs_thresh=cs_thresh,zmagcut=zmagcut,
         imaflags_thresh=imaflags_thresh)
   num = len(zcat)
   #indexes = compress(crit,arange(len(crit)))

   # write to dropout catalog
   # process headers
   HEADER = ""
   offset = 0
   HEADER += "# 1 ID_MOSAIC\n"
   HEADER += "# 2 ALPHA_J2000\n"
   HEADER += "# 3 DELTA_J2000\n"
   HEADER += "# 4 SECT_REFNUM\n"
   HEADER += "# 5 THETA_IMAGE\n"
   HEADER += "# 6 ELLIPTICITY\n"
   HEADER += "# 7 X_SECT\n"
   HEADER += "# 8 Y_SECT\n"
   HEADER += "# 9 B_FLUX_RADIUS_1\n"
   HEADER += "# 10 V_FLUX_RADIUS_1\n"
   HEADER += "# 11 I_FLUX_RADIUS_1\n"
   HEADER += "# 12 Z_FLUX_RADIUS_1\n"
   HEADER += "# 13 B_MAG_AUTO\n"
   HEADER += "# 14 V_MAG_AUTO\n"
   HEADER += "# 15 I_MAG_AUTO\n"
   HEADER += "# 16 Z_MAG_AUTO\n"
   HEADER += "# 17 PHOTZ_WEIGHTED\n"
   HEADER += "# 18 PHOTZ_95MIN\n"
   HEADER += "# 19 PHOTZ_95MAX\n"
   HEADER += "# 20 SPECZ\n"
   HEADER += "# 21 CULLFLAG\n"
   HEADER += "# 22 B_MAG_ISO\n"
   HEADER += "# 23 B_MAGERR_ISO\n"
   HEADER += "# 24 V_MAG_ISO\n"
   HEADER += "# 25 V_MAGERR_ISO\n"
   HEADER += "# 26 I_MAG_ISO\n"
   HEADER += "# 27 I_MAGERR_ISO\n"
   HEADER += "# 28 Z_MAG_ISO\n"
   HEADER += "# 29 Z_MAGERR_ISO\n"
   HEADER += "# 30 B_FLUX_ISO\n"
   HEADER += "# 31 B_FLUXERR_ISO\n"
   HEADER += "# 32 V_FLUX_ISO\n"
   HEADER += "# 33 V_FLUXERR_ISO\n"
   HEADER += "# 34 I_FLUX_ISO\n"
   HEADER += "# 35 I_FLUXERR_ISO\n"
   HEADER += "# 36 Z_FLUX_ISO\n"
   HEADER += "# 37 Z_FLUXERR_ISO\n"
   HEADER += "# 38 B_FLUX_AUTO\n"
   HEADER += "# 39 B_FLUXERR_AUTO\n"
   HEADER += "# 40 V_FLUX_AUTO\n"
   HEADER += "# 41 V_FLUXERR_AUTO\n"
   HEADER += "# 42 I_FLUX_AUTO\n"
   HEADER += "# 43 I_FLUXERR_AUTO\n"
   HEADER += "# 44 Z_FLUX_AUTO\n"
   HEADER += "# 45 Z_FLUXERR_AUTO\n"
   HEADER += "# 46 B_IMAFLAGS_ISO\n"
   HEADER += "# 47 V_IMAFLAGS_ISO\n"
   HEADER += "# 48 I_IMAFLAGS_ISO\n"
   HEADER += "# 49 Z_IMAFLAGS_ISO\n"
   HEADER += "# 50 Z_CLASS_STAR\n"


   f = open(outname,'w')
   #f.write(HEADER)
   # construct arrays with selected objects only
   id = compress(crit,zcat.id_mosaic)
   ra = compress(crit,zcat.alpha_j2000)
   dec = compress(crit,zcat.delta_j2000)
   sect = compress(crit,zcat.sect_refnum)
   theta = compress(crit,zcat.theta_image)
   ell = compress(crit,zcat.ellipticity)
   ximage = compress(crit,zcat.x_sect)
   yimage = compress(crit,zcat.y_sect)
   bflux_radius_1 = compress(crit,bcat.flux_radius_1)
   vflux_radius_1 = compress(crit,vcat.flux_radius_1)
   iflux_radius_1 = compress(crit,icat.flux_radius_1)
   zflux_radius_1 = compress(crit,zcat.flux_radius_1)
   bmag_auto = compress(crit,bcat.mag_auto)
   vmag_auto = compress(crit,vcat.mag_auto)
   imag_auto = compress(crit,icat.mag_auto)
   zmag_auto = compress(crit,zcat.mag_auto)
   bmag_iso = compress(crit, bcat.mag_iso)
   bmagerr_iso = compress(crit, bcat.magerr_iso)
   vmag_iso = compress(crit, vcat.mag_iso)
   vmagerr_iso = compress(crit, vcat.magerr_iso)
   imag_iso = compress(crit, icat.mag_iso)
   imagerr_iso = compress(crit, icat.magerr_iso)
   zmag_iso = compress(crit, zcat.mag_iso)
   zmagerr_iso = compress(crit, zcat.magerr_iso)
   bflux_iso = compress(crit, bcat.flux_iso)
   bfluxerr_iso = compress(crit, bcat.fluxerr_iso)
   vflux_iso = compress(crit, vcat.flux_iso)
   vfluxerr_iso = compress(crit, vcat.fluxerr_iso)
   iflux_iso = compress(crit, icat.flux_iso)
   ifluxerr_iso = compress(crit, icat.fluxerr_iso)
   zflux_iso = compress(crit, zcat.flux_iso)
   zfluxerr_iso = compress(crit, zcat.fluxerr_iso)
   bflux_auto = compress(crit, bcat.flux_auto)
   bfluxerr_auto = compress(crit, bcat.fluxerr_auto)
   vflux_auto = compress(crit, vcat.flux_auto)
   vfluxerr_auto = compress(crit, vcat.fluxerr_auto)
   iflux_auto = compress(crit, icat.flux_auto)
   ifluxerr_auto = compress(crit, icat.fluxerr_auto)
   zflux_auto = compress(crit, zcat.flux_auto)
   zfluxerr_auto = compress(crit, zcat.fluxerr_auto)
   bflag = compress(crit, bcat.imaflags_iso)
   vflag = compress(crit, vcat.imaflags_iso)
   iflag = compress(crit, icat.imaflags_iso)
   zflag = compress(crit, zcat.imaflags_iso)
   class_star = compress(crit, zcat.class_star)


   # match with photz/specz catalog in GOODS-S
   if max(dec) < 0.:  # if in South
      cz = sextractor(photzcat)
      photz = zeros(len(id))
      photz_95min = zeros(len(id))
      photz_95max = zeros(len(id))
      specz = zeros(len(id))
      for i in range(len(id)):
         if id[i] in cz.id_zband:
            photz[i] = compress(cz.id_zband==id[i],cz.photz_weighted)[0]
            photz_95min[i] = compress(cz.id_zband==id[i],cz.photz_95min)[0]
            photz_95max[i] = compress(cz.id_zband==id[i],cz.photz_95max)[0]
            specz[i] = compress(cz.id_zband==id[i],cz.specz)[0]
         else:
            photz[i]=-9.0; photz_95min[i]=-9.0; photz_95max[i]=-9.0; specz[i]=-9.0
   else:
      photz = ones(len(id))*(-9.0)
      photz_95min = ones(len(id))*(-9.0)
      photz_95max = ones(len(id))*(-9.0)
      specz = ones(len(id))*(-9.0)

   # match with spec-z catalog in GOODS-N
   if min(dec)>0.: # if in GOODS-N
      n_specz = sextractor('/Users/khuang/Dropbox/Research/catalogs/specz_photz/goods_n_specz_0111.txt')
      tol_deg = tol_as/3600.
      for i in range(len(id)):
         angdist = angsep.angsep(n_specz.ra,n_specz.dec,ra[i],dec[i])
         if min(angdist) <= tol_deg:
            j = argsort(angdist)[0]
            specz[i] = n_specz.specz[j]

   # match with existing cull flag file
   cc = sextractor(cullcat)
   cullflag = ones(len(id),'int')*(-1)
   for i in range(len(id)):
      for j in range(len(cc)):
         if (id[i]==cc.id_mosaic[j]) & (dec[i]*cc.delta_j2000[j]>0.):
            cullflag[i] = cc.cullflag[j]
   
   f.write(HEADER)
   # write entries
   for i in range(len(id)):
      f.write('%d ' % id[i])  # ACS ID
      f.write('%f ' % ra[i])  # RA
      f.write('%f ' % dec[i])  # DEC
      f.write('%d ' % sect[i])  # SECT NUMBER
      f.write('%f ' % theta[i])  # THETA_IMAGE
      f.write('%f ' % ell[i])  # ELLIPTICITY
      f.write('%f ' % ximage[i])  # X_SECT
      f.write('%f ' % yimage[i])  # Y_SECT
      f.write('%f ' % bflux_radius_1[i]) # half-light radius in SExtractor
      f.write('%f ' % vflux_radius_1[i])
      f.write('%f ' % iflux_radius_1[i])
      f.write('%f ' % zflux_radius_1[i])
      f.write('%f ' % bmag_auto[i]) # mag_auto
      f.write('%f ' % vmag_auto[i])
      f.write('%f ' % imag_auto[i])
      f.write('%f ' % zmag_auto[i])
      f.write('%f ' % photz[i])
      f.write('%f ' % photz_95min[i])
      f.write('%f ' % photz_95max[i])
      f.write('%f ' % specz[i])
      f.write('%d ' % cullflag[i]) 
      f.write('%f %f ' % (bmag_iso[i], bmagerr_iso[i])) # 22-23
      f.write('%f %f ' % (vmag_iso[i], vmagerr_iso[i])) # 24-25
      f.write('%f %f ' % (imag_iso[i], imagerr_iso[i])) # 26-27
      f.write('%f %f ' % (zmag_iso[i], zmagerr_iso[i])) # 28-29
      f.write('%f %f ' % (bflux_iso[i], bfluxerr_iso[i])) # 30-31
      f.write('%f %f ' % (vflux_iso[i], vfluxerr_iso[i])) # 32-33
      f.write('%f %f ' % (iflux_iso[i], ifluxerr_iso[i])) # 34-35
      f.write('%f %f ' % (zflux_iso[i], zfluxerr_iso[i])) # 36-37
      f.write('%f %f ' % (bflux_auto[i], bfluxerr_auto[i])) # 38-39
      f.write('%f %f ' % (vflux_auto[i], vfluxerr_auto[i])) # 40-41
      f.write('%f %f ' % (iflux_auto[i], ifluxerr_auto[i])) # 42-43
      f.write('%f %f ' % (zflux_auto[i], zfluxerr_auto[i])) # 44-45
      f.write('%d %d ' % (bflag[i], vflag[i])) # 46-47
      f.write('%d %d ' % (iflag[i], zflag[i])) # 48-49
      f.write('%f ' % class_star[i]) # 50
      f.write('\n')
   f.flush()
   f.close()
   print "B-dropout catalog made; total # =", len(id)

def make_bdropscat_all(northpar,southpar,outcat):
   # Make B-dropout catalogs for both North and South, then combine them

   # North first
   cn = {}
   cn = parseconfig(northpar,paramdict=cn)
   bcat,vcat,icat,zcat=cn['BCAT'],cn['VCAT'],cn['ICAT'],cn['ZCAT']
   outname_n=cn['OUTNAME']
   photzcat=None
   cullcat=cn['CULLCAT']
   cs_thresh=cn['CS_THRESH']
   zmagcut = cn['ZMAGCUT']
   imaflags_thresh = cn['IMAFLAGS_THRESH']
   try:
      magform = cn['MAGFORM']
   except: magform = 'magiso'

   make_bdropcat(bcat,vcat,icat,zcat,outname_n,photzcat,cullcat,cs_thresh=cs_thresh,
      zmagcut=zmagcut,imaflags_thresh=imaflags_thresh,magform=magform,tol_as=0.5)

   # South second
   cs = {}
   cs = parseconfig(southpar,paramdict=cs)
   bcat,vcat,icat,zcat=cs['BCAT'],cs['VCAT'],cs['ICAT'],cs['ZCAT']
   outname_s=cs['OUTNAME']
   photzcat=cs['PHOTZCAT']
   cullcat=cs['CULLCAT']
   cs_thresh=cs['CS_THRESH']
   zmagcut = cs['ZMAGCUT']
   imaflags_thresh = cs['IMAFLAGS_THRESH']
   try:
      magform = cs['MAGFORM']
   except: magform = 'magiso'

   make_bdropcat(bcat,vcat,icat,zcat,outname_s,photzcat,cullcat,cs_thresh=cs_thresh,
      zmagcut=zmagcut,imaflags_thresh=imaflags_thresh,magform=magform,tol_as=0.5)

   # Now merge north & south
   f = open(outcat,"w")
   fn = open(outname_n)
   fs = open(outname_s)
   lines_n = fn.readlines()
   for i in range(len(lines_n)):
      f.write(lines_n[i])
   lines_s = fs.readlines()
   for i in range(len(lines_s)):
      if lines_s[i].startswith("#"): pass
      else: f.write(lines_s[i])
   f.close()
   return 0

# for B-dropout selection, need B-, V-, and z-band catalogs
bcat_n = '/Users/khuang/Dropbox/Research/catalogs/goods/north/acsz_v2/h_goods_nb_r2.0z.cat'
vcat_n = '/Users/khuang/Dropbox/Research/catalogs/goods/north/acsz_v2/h_goods_nv_r2.0z.cat'
icat_n = '/Users/khuang/Dropbox/Research/catalogs/goods/north/acsz_v2/h_goods_ni_r2.0z.cat'
zcat_n = '/Users/khuang/Dropbox/Research/catalogs/goods/north/acsz_v2/h_goods_nz_r2.0z.cat'
bcat_s = '/Users/khuang/Dropbox/Research/catalogs/goods/south/acsz_v2/h_goods_sb_r2.0z.cat'
vcat_s = '/Users/khuang/Dropbox/Research/catalogs/goods/south/acsz_v2/h_goods_sv_r2.0z.cat'
icat_s = '/Users/khuang/Dropbox/Research/catalogs/goods/south/acsz_v2/h_goods_si_r2.0z.cat'
zcat_s = '/Users/khuang/Dropbox/Research/catalogs/goods/south/acsz_v2/h_goods_sz_r2.0z.cat'


def b_readcats(bcatloc,vcatloc,zcatloc):
   bcat = sextractor(bcatloc)
   vcat = sextractor(vcatloc)
   zcat = sextractor(zcatloc)
   return bcat,vcat,zcat


if __name__ == "__main__":
   print "Make B-dropout catalog for GOODS North or South"
   parfile = sys.argv[1]
   c = parseconfig(parfile)
   bcat = c['BCAT']
   vcat = c['VCAT']
   icat = c['ICAT']
   zcat = c['ZCAT']
   outname = c['OUTNAME']
   try:
      photzcat = c['PHOTZCAT']
   except:
      photzcat = None
   cs_thresh = c['CS_THRESH']
   zmagcut = c['ZMAGCUT']
   imaflags_thresh = c['IMAFLAGS_THRESH']
   try:
      magform = c['MAGFORM']
   except: magform = 'magiso'
   cullcat = c['CULLCAT']
   make_bdropcat(bcat,vcat,icat,zcat,outname,photzcat,cullcat,cs_thresh=cs_thresh,
      zmagcut=zmagcut,imaflags_thresh=imaflags_thresh,magform=magform)
   
