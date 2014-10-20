#!/usr/bin/env python

from numpy import *
from pygoods import *
import os, sys, glob



#def colorsel_g04(vmags,imags,zmags):
#   """Giavalisco et al. 2004 V-dropouts color selection
#      colorsel_g04(vmags,imags,zmags)
#   """
#   vmi = vmags - imags
#   imz = imags - zmags
#   c1 = (vmi > 1.5 + 0.9 * imz)
#   c2 = (vmi > 2.0)
#   c3 = (vmi >= 1.2)
#   c4 = (imz <= 1.3)
#   return  ((c1+c2) * c3 * c4)
def colorsel_g04(vflux,vfluxerr,vzpt,iflux,ifluxerr,izpt,zflux,zfluxerr,zzpt):
	"""Giavalisco et al. 2004 criteria"""
	c0 = (zflux / zfluxerr >= 5.0)  # z-band S/N criteria
	vmags = flux2mag_ul(vflux,vfluxerr,vzpt)
	imags = flux2mag_ul(iflux,ifluxerr,izpt)
	zmags = flux2mag(zflux,zfluxerr,zzpt)
	vmi = vmags - imags
	imz = imags - zmags
	c1 = (vmi > 1.5 + 0.9 * imz)
	c2 = (vmi > 2.0)
	c3 = (vmi >= 1.2)
	c4 = (imz <= 1.3)
	return ((c1 | c2) & c3 & c4)


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

def dropsel_magauto(bcat,vcat,icat,zcat,cs_thresh=0.95,zmagcut=26.5,
   imaflags_thresh=4,bstonlim=2.0):
   """Select V-dropouts with following requirements:
      - satisfy color criteria
      - class_star <= cs_thresh in i, z
      - imaflags_iso <= imaflags_thresh in V, i, z
      - zmag <= zmagcut
      - S/N <= 2.0 in B-band
      - S/N >= 5.0 in z-band
   """
   # class_star threshold
   cs_crit = (zcat.class_star<=cs_thresh)+(zcat.mag_auto>26.0)
   # z-band magnitude cut
   zmag_crit = (zcat.mag_auto<=zmagcut)
   # image flags
   flag_crit = (vcat.imaflags_iso<4)*(icat.imaflags_iso<4)*(zcat.imaflags_iso<4)
   # S/N cut - clarification on how to calculate S/N?
   zston_crit = (abs(zcat.flux_auto/zcat.fluxerr_auto)>=5.)
   bston_crit = (abs(bcat.flux_auto/bcat.fluxerr_auto)<=bstonlim)

   # color selection
   vmags = vcat.mag_auto.copy()
   imags = icat.mag_auto.copy()
   zmags = zcat.mag_auto.copy()
   # Use 1-sigma as upper limit in V-band if S/N <= 1.0
   for i in range(len(vcat)):
      if vcat.flux_auto[i]/vcat.fluxerr_auto[i] <= 1.0:
         vmags[i] = flux2mag(vcat.fluxerr_auto[i],26.49341)
      if icat.flux_auto[i]/icat.fluxerr_auto[i] <= 1.0:
         imags[i] = flux2mag(icat.fluxerr_auto[i],25.64053)
      if zcat.flux_auto[i]/zcat.fluxerr_auto[i] <= 1.0:
         zmags[i] = flux2mag(zcat.fluxerr_auto[i],24.84315)

   color_crit = colorsel_g04(vmags,imags,zmags)

   crit = cs_crit*zmag_crit*flag_crit*zston_crit*bston_crit*color_crit
   return crit


def dropsel_magiso(bcat,vcat,icat,zcat,cs_thresh=0.95,zmagcut=26.5,
   imaflags_thresh=4,bstonlim=2.0):
   """Select V-dropouts with following requirements:
      - satisfy color criteria
      - class_star <= cs_thresh in i, z
      - imaflags_iso <= imaflags_thresh in V, i, z
      - zmag <= zmagcut
      - S/N <= 2.0 in B-band
      - S/N >= 5.0 in z-band
   """
   # class_star threshold
   cs_crit = (zcat.class_star<=cs_thresh)+(zcat.mag_auto>26.0)
   # z-band magnitude cut
   zmag_crit = (zcat.mag_auto<=zmagcut)
   # image flags
   flag_crit = (icat.imaflags_iso<imaflags_thresh)&(zcat.imaflags_iso<imaflags_thresh)
   # S/N cut - clarification on how to calculate S/N?
   zston_crit = (abs(zcat.flux_auto/zcat.fluxerr_auto)>=5.)
   bston_crit = ((bcat.flux_auto/bcat.fluxerr_auto)<=bstonlim)

   # color selection
   vmags = vcat.mag_iso.copy()
   imags = icat.mag_iso.copy()
   zmags = zcat.mag_iso.copy()
   # Use 1-sigma as upper limit in V-band if S/N <= 1.0
   for i in range(len(vcat)):
      if vcat.flux_iso[i]/vcat.fluxerr_iso[i] <= 1.0:
         vmags[i] = flux2mag(vcat.fluxerr_iso[i],26.49341)
      if icat.flux_iso[i]/icat.fluxerr_iso[i] <= 1.0:
         imags[i] = flux2mag(icat.fluxerr_iso[i],25.64053)
      if zcat.flux_iso[i]/zcat.fluxerr_iso[i] <= 1.0:
         zmags[i] = flux2mag(zcat.fluxerr_iso[i],24.84315)

   color_crit = colorsel_g04(vmags,imags,zmags)

   crit = cs_crit & zmag_crit & color_crit & zston_crit & bston_crit & flag_crit
   return crit


def make_vdropcat(bcat,vcat,icat,zcat,outname,photzcat,cullcat,cs_thresh=0.95,zmagcut=26.5,
   imaflags_thresh=4,magform='magiso',tol_as=0.5,bstonlim=2.0):
   # makes V-dropout catalog in either North or South
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
      crit = dropsel_magauto(bcat,vcat,icat,zcat,cs_thresh=0.95,zmagcut=zmagcut,
         imaflags_thresh=imaflags_thresh)
   elif magform == 'magiso':
      crit = dropsel_magiso(bcat,vcat,icat,zcat,cs_thresh=0.95,zmagcut=zmagcut,
         imaflags_thresh=imaflags_thresh,bstonlim=bstonlim)
   num = len(zcat)
   #indexes = compress(crit,arange(len(crit)))

   # write to dropout catalog
   # process headers
   HEADER = ""
   offset = 0
   HEADER += "# 1 ID_MOSAIC\n"
   HEADER += "# 2 ALPHA_J2000\n"
   HEADER += "# 3 DELTA_J2000\n"
   HEADER += "# 4 OSECT_REFNUM\n"
   HEADER += "# 5 THETA_IMAGE\n"
   HEADER += "# 6 ELLIPTICITY\n"
   HEADER += "# 7 X_IMAGE\n"
   HEADER += "# 8 Y_IMAGE\n"
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
   id_drop = compress(crit,zcat.id_mosaic)
   ra = compress(crit,zcat.alpha_j2000)
   dec = compress(crit,zcat.delta_j2000)
   osect = compress(crit,zcat.osect_refnum)
   theta = compress(crit,zcat.theta_image)
   ell = compress(crit,zcat.ellipticity)
   ximage = compress(crit,zcat.x_image)
   yimage = compress(crit,zcat.y_image)
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
         else: photz[i]=-9.0; photz_95min[i]=-9.0; photz_95max[i]=-9.0; specz[i]=-9.0
   else:
      photz = ones(len(id))*(-9.0)
      photz_95min = ones(len(id))*(-9.0)
      photz_95max = ones(len(id))*(-9.0)
      specz = ones(len(id))*(-9.0)

   # match with spec-z catalog in GOODS-N
   if min(dec)>0.: # if in GOODS-N
      n_photz = sextractor('/Users/kuang/Dropbox/Research/catalogs/specz_photz/z.hdfn.r2.0z.20090224_sex.dat')
      tol_deg = tol_as/3600.
      for i in range(len(id)):
         angdist = angsep.angsep(n_photz.ra,n_photz.dec,ra[i],dec[i])
         if min(angdist) <= tol_deg:
            j = argsort(angdist)[0]
            specz[i] = n_photz.specz[j]

   # match with existing cull flag file
   cc = sextractor(cullcat)
   cullflag = ones(len(id),'int')*(-1)
   for i in range(len(id)):
      for j in range(len(cc)):
         if (id[i]==cc._1[j])*(dec[i]*cc._3[j]>0.):
            cullflag[i] = cc._4[j]

   f.write(HEADER)
   # write entries
   for i in range(len(id)):
      f.write('%d ' % id[i])  # 1
      f.write('%f ' % ra[i])  # 2
      f.write('%f ' % dec[i])  # 3
      f.write('%d ' % osect[i])  # 4
      f.write('%f ' % theta[i])  # 5
      f.write('%f ' % ell[i])  # 6
      f.write('%f ' % ximage[i])  # 7
      f.write('%f ' % yimage[i])  # 8
      f.write('%f ' % bflux_radius_1[i]) # 9
      f.write('%f ' % vflux_radius_1[i]) # 10
      f.write('%f ' % iflux_radius_1[i]) # 11
      f.write('%f ' % zflux_radius_1[i]) # 12
      f.write('%f ' % bmag_auto[i])  # 13
      f.write('%f ' % vmag_auto[i])  # 14
      f.write('%f ' % imag_auto[i])  # 15
      f.write('%f ' % zmag_auto[i])  # 16
      f.write('%f ' % photz[i])  # 17
      f.write('%f ' % photz_95min[i])  # 18
      f.write('%f ' % photz_95max[i])  # 19
      f.write('%f ' % specz[i])  # 20
      f.write('%d ' % cullflag[i])  # 21
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
   print "V-dropout catalog made; total # =", len(id)


def make_vdropscat_all(northpar,southpar,outcat):
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
   bn = cn['BSTONLIM']
   try:
      magform = cn['MAGFORM']
   except: magform = 'magiso'

   make_vdropcat(bcat,vcat,icat,zcat,outname_n,photzcat,cullcat,cs_thresh=cs_thresh,
      zmagcut=zmagcut,imaflags_thresh=imaflags_thresh,magform=magform,tol_as=0.5,
      bstonlim=bn) 

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
   bs = cs['BSTONLIM']
   try:
      magform = cs['MAGFORM']
   except: magform = 'magiso'

   make_vdropcat(bcat,vcat,icat,zcat,outname_s,photzcat,cullcat,cs_thresh=cs_thresh,
      zmagcut=zmagcut,imaflags_thresh=imaflags_thresh,magform=magform,tol_as=0.5,
      bstonlim=bs)

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

# for V-dropout selection, need B-, V-, i-, and z-band catalogs
bcat_n = '/data/raid4/catalogs/team/n/acsz_v1.9/photom/nb_acs_d1.0z.cat'
vcat_n = '/data/raid4/catalogs/team/n/acsz_v1.9/photom/nv_acs_d1.0z.cat'
icat_n = '/data/raid4/catalogs/team/n/acsz_v1.9/photom/ni_acs_d1.0z.cat'
zcat_n = '/data/raid4/catalogs/team/n/acsz_v1.9/photom/nz_acs_d1.0z.cat'
bcat_s = '/data/raid4/catalogs/team/s/acsz_v1.9/photom/sb_acs_d1.0z.cat'
vcat_s = '/data/raid4/catalogs/team/s/acsz_v1.9/photom/sv_acs_d1.0z.cat'
icat_s = '/data/raid4/catalogs/team/s/acsz_v1.9/photom/si_acs_d1.0z.cat'
zcat_s = '/data/raid4/catalogs/team/s/acsz_v1.9/photom/sz_acs_d1.0z.cat'

dc_header = []
dc_header += ["FLUX_RADIUS"]
dc_header += ["FLUX_RADIUS_1"]
dc_header += ["FLUX_RADIUS_2"]
dc_header += ["CLASS_STAR"]
dc_header += ["MAG_ISO"]
dc_header += ["MAGERR_ISO"]
dc_header += ["MAG_AUTO"]
dc_header += ["MAGERR_AUTO"]

def readcats(vcatloc,icatloc,zcatloc):
   bcat = sextractor(bcatloc)
   vcat = sextractor(vcatloc)
   icat = sextractor(icatloc)
   zcat = sextractor(zcatloc)
   return bcat, vcat, icat, zcat

if __name__ == "__main__":
   print "Make V-dropout catalog for GOODS North or South"
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
   cullcat = c['CULLCAT']
   make_vdropcat(bcat,vcat,icat,zcat,outname,photzcat,cullcat,cs_thresh=cs_thresh,\
      zmagcut=zmagcut,imaflags_thresh=imaflags_thresh)
