#!/usr/bin/env python

import numpy as np
import os, sys, glob
import subprocess
from pygoods import sextractor, Ftable, fitstable
import fitsutil
import pyfits
import yaml
from sedtools.bands import filters
from datetime import datetime
from astropy import coordinates

sex_exec = '/Users/khuang/bin/cex'
surfsup_dir = '/Users/khuang/Dropbox/codes/mypytools/surfsup'
clash_bands = ['f435w','f475w','f606w','f625w','f775w','f814w','f850lp','f105w','f110w','f125w','f140w','f160w']

def magerr_from_sn(sn):
  return 2.5 * np.log10(1.+1./sn)

def sn_from_magerr(magerr):
  return 1. / (10. ** (0.4 * magerr) - 1.)

def abmag_to_uJy(mag):
  return 10. ** (-0.4 * (mag - 23.9))

class HSTPipeline(object):
  def __init__(self, cluster_name, paramfile_name):
    paramfile = os.path.join(surfsup_dir, paramfile_name)
    if not os.path.exists(paramfile):
      print "Parameter file %s does not exists." % os.path.split(paramfile)[-1]
      print "The parameter files are available for the following clusters:"
      pfiles = glob.glob(surfsup_dir+'/*_hst_photom.dat.yml')
      for p in pfiles:
        clname = os.path.split(p)[-1].split('_')[0].upper()
        print clname
      raise ValueError
    c = yaml.load(open(paramfile, 'rb'))
    self.c = c
    scale = c['scale'] * 1000  # in milli-arcsec
    scale = int(round(scale))
    scalestr = '%2dmas' % scale
    for k in c.keys():
      setattr(self, k, c[k])
    self.cluster_name = cluster_name
    self.fixcols = ['number', 'x_image', 'y_image', 'alpha_j2000', 'delta_j2000']
    for f in self.filters:
      self.zeropoints[f] = self.zeropoints[f] - self.ebmv_mean * self.ext_coeff[f]
    self.sex_config = {}
    for f in self.filters:
      self.sex_config[f] = '%s_%s_%s.sex' % (f, self.cluster_name.lower(), scalestr)
      os.system('cp %s %s' % (self.sexfile, self.sex_config[f]))
    self.sex_catalog = {}
    for f in self.filters:
      self.sex_catalog[f] = '%s_%s_%s.cat' % (f, self.cluster_name.lower(), scalestr)
    self.sex_fitstable = {}
    for f in self.filters:
      self.sex_fitstable[f] = '%s_%s_%s.fits' % (f, self.cluster_name.lower(), scalestr)
      if os.path.exists(self.sex_fitstable[f].lower()):
        os.remove(self.sex_fitstable[f])
    try:
      self.wht_type = c['wht_type'].upper()
    except:
      self.wht_type = 'MAP_WEIGHT'

  def make_flag(self, overwrite=False):
    """
    Make flag images from input weight images. Assume that input weight images
    have file names like *_wht.fits.
    """
    for f in self.filters:
      w = self.wht_images[f]
      hdr = pyfits.getheader(w)
      wht_array = pyfits.getdata(w)
      flg_name = w.replace('wht', 'flg')
      if (overwrite==True) or (os.path.exists(flg_name)==False):
        flg_array = np.where(wht_array > 0, 0, 1).astype('int16')
        pyfits.append(flg_name, flg_array, hdr)
        print "Made flag image from %s..." % w

  def run_sextractor(self, filters=[], **custom_args):
    # Run SExtractor for each band
    if not len(filters):
      filters = self.filters
    for f in filters:
      args = ['%s,%s' % (self.drz_images[self.detectband], self.drz_images[f])]
      args += ['-c', self.sex_config[f]]
      args += ['-CATALOG_NAME', self.sex_catalog[f]]
      if f == self.detectband:
        args += ['-CHECKIMAGE_NAME', self.drz_images[f].replace('drz', 'seg')]
      args += ['-BACKPHOTO_TYPE', 'LOCAL']
      args += ['-WEIGHT_IMAGE', '%s,%s' % (self.wht_images[self.detectband], self.wht_images[f])]
      args += ['-WEIGHT_TYPE', '%s,%s' % (self.wht_type, self.wht_type)]
      args += ['-FLAG_IMAGE', self.flg_images[f]]
      args += ['-MAG_ZEROPOINT', str(self.zeropoints[f])]
      args += ['-GAIN', str(self.effgain[f])]
      for cargs in custom_args.keys():
        args += ['-%s' % cargs.upper(), str(custom_args[cargs])]
      args = [sex_exec] + args
      # print ' '.join(args)
      print args
      subprocess.call(args)
      self.sex2fits(f)
      self.add_filternames(f)

  def sex2fits(self, f, overwrite=True):
    # Convert SExtractor catalogs into FITS tables
    print "Convert SExtractor catalogs into FITS tables..."
    if (overwrite==True) or (os.path.exists(self.sex_fitstable[f]==False)):
      c = sextractor(self.sex_catalog[f])
      fitsutil.sex2fits(c, self.sex_fitstable[f])

  def add_filternames(self, f):
    # prepend the FITS table columns by filter name
    print "Adding filter names to the columns..."
    c = Ftable(self.sex_fitstable[f])
    old_colnames = []
    new_colnames = []
    for col in c.Columns:
      if col.lower() in self.fixcols:
        pass
      else:
        old_colnames += [col]
        new_colnames += ['%s_%s' % (f, col)]
    fitsutil.change_column_names(self.sex_fitstable[f], old_colnames, 
                               new_colnames)

  def merge_sexcats(self):
    """
    Assume that files self.sex_catalog exist.
    """
    # Match & merge photometry into one master catalog
    # Note: update all magnitude errors with the scaled flux errors
    merged_columns = []
    for f in self.filters:
      c = Ftable(self.sex_fitstable[f])
      skyrms_factor = self.c['skyrms_factor'][f]
      if f == self.detectband:
        merged_columns += [pyfits.Column(name='number', array=c.number, format='I')]
        merged_columns += [pyfits.Column(name='x_image', array=c.x_image, format='D')]
        merged_columns += [pyfits.Column(name='y_image', array=c.y_image, format='D')]
        merged_columns += [pyfits.Column(name='alpha_j2000', array=c.alpha_j2000, format='D')]
        merged_columns += [pyfits.Column(name='delta_j2000', array=c.delta_j2000, format='D')]
      flux_iso = c.__getitem__('%s_flux_iso' % f)
      fluxerr_iso = c.__getitem__('%s_fluxerr_iso' % f)
      print "Scaling fluxerr_iso..."
      fluxerr_iso_scaled = fluxerr_iso * skyrms_factor
      # sn_iso = flux_iso / fluxerr_iso
      sn_iso = flux_iso / fluxerr_iso_scaled
      flux_auto = c.__getitem__('%s_flux_auto' % f)
      fluxerr_auto = c.__getitem__('%s_fluxerr_auto' % f)
      sn_auto = flux_auto / fluxerr_auto
      print "Scaling fluxerr_auto..."
      fluxerr_auto_scaled = fluxerr_auto * skyrms_factor
      for i in range(len(c.Columns)):
        colname = c.Columns[i]
        if colname in self.fixcols:
          continue
        else:
          merged_columns += [pyfits.Column(name=colname, 
                            array=c.__getitem__(colname),
                            format=c.d.formats[i])]
          if colname == '%s_fluxerr_iso' % f:
            # also add the scaled fluxerr 
            merged_columns += [pyfits.Column(name='%s_fluxerr_iso_scaled' % f,
                               array=fluxerr_iso_scaled,
                               format='D')]
            magerr_iso_scaled = magerr_from_sn(sn_iso)
          if colname == '%s_fluxerr_auto' % f:
            merged_columns += [pyfits.Column(name='%s_fluxerr_auto_scaled' % f,
                                array=fluxerr_auto_scaled,
                                format='D')]
            magerr_auto_scaled = magerr_from_sn(sn_auto)
          if colname.startswith('%s_fluxerr_aper' % f):
            fluxerr_aper_n = getattr(c, colname)
            print "Scaling %s..." % colname
            fluxerr_aper_n_scaled = fluxerr_aper_n * skyrms_factor
            flux_aper_n = getattr(c, '%s_flux_aper' % f)
            sn_aper = flux_aper_n / fluxerr_aper_n_scaled
            merged_columns += [pyfits.Column(name='%s_scaled' % colname,
                               array=fluxerr_aper_n_scaled,
                               format='D')]
          # Now update MAG_ISO: MAG_ISO=99.0 if S/N < 1 and MAGERR_ISO becomes 
          # the 1-sigma magnitude limit; MAG_ISO=-99.0 if object is not detected
          if colname == '%s_mag_iso' % f:
            mag_iso = np.where(sn_iso >= 1., c.__getitem__('%s_mag_iso' % f),
                               99.0)
            mag_iso = np.where(fluxerr_iso==0, -99.0, mag_iso)
            merged_columns[-1].array = mag_iso   # update the array
          if colname == '%s_magerr_iso' % f:
            magerr_iso = np.where(fluxerr_iso_scaled==0, 0., 
                                  2.5*np.log10(1.+1./sn_iso))
            # if S/N <= 1., use the 1-sigma magnitude limit as magerr_iso
            magerr_iso = np.where(sn_iso>=1., magerr_iso, 
                          self.zeropoints[f]-2.5*np.log10(fluxerr_iso_scaled))
            merged_columns[-1].array = magerr_iso
          # Write the 1-sigma MAG_ISO for ease of color selection
          if colname == '%s_mag_auto' % f:
            mag_auto = np.where(sn_auto >= 1., c.__getitem__('%s_mag_auto' % f),
                                99.0)
            mag_auto = np.where(fluxerr_auto==0., -99.0, mag_auto)
            merged_columns[-1].array = mag_auto
          if colname == '%s_magerr_auto' % f:
            magerr_auto = np.where(fluxerr_auto_scaled==0, 0., 
                                   2.5*np.log10(1.+1./sn_auto))
            magerr_auto = np.where(sn_auto>=1., magerr_auto,
                          self.zeropoints[f]-2.5*np.log10(fluxerr_auto_scaled))
      mag_iso_1sig = np.where(sn_iso>=1.0, mag_iso, magerr_iso)
      merged_columns += [pyfits.Column(name='%s_mag_iso_1sig' % f,
                         array=mag_iso_1sig, format='D')]
      mag_auto_1sig = np.where(sn_auto>=1.0, mag_auto, magerr_auto)
      merged_columns += [pyfits.Column(name='%s_mag_auto_1sig' % f,
                          array=mag_auto_1sig, format='D')]
      naper = 0
      flux_aper = c.__getitem__('%s_flux_aper' % f)
      fluxerr_aper = c.__getitem__('%s_fluxerr_aper' % f)
      print "Scaling fluxerr_aper..."
      fluxerr_aper_scaled = fluxerr_aper * skyrms_factor
      # sn_aper = flux_aper / fluxerr_aper
      sn_aper = flux_aper / fluxerr_aper_scaled
      for i in range(len(merged_columns)):
        if merged_columns[i].name == '%s_mag_aper' % f:
          mag_aper = np.where(sn_aper >= 1., c.__getitem__('%s_mag_aper'%f),
                              99.0)
          mag_aper = np.where(fluxerr_aper==0, -99.0, mag_aper)
          merged_columns[i].array = mag_aper
        if merged_columns[i].name == '%s_magerr_aper' % f:
          # magerr_aper = np.where(fluxerr_aper==0, 0., 
          #                        c.__getitem__('%s_magerr_aper'%f))
          # magerr_aper = np.where(sn_aper>=1., magerr_aper, 
          #                        self.zeropoints[f]-2.5*np.log10(fluxerr_aper))
          magerr_aper = magerr_from_sn(sn_aper)
          magerr_aper = np.where(fluxerr_aper_scaled==0., 0.,
                                 2.5*np.log10(1.+1./sn_aper))
          # update magerr_aper with 1-sigma magnitude limit
          magerr_aper = np.where(sn_aper>=1., magerr_aper, 
                                 self.zeropoints[f]-2.5*np.log10(fluxerr_aper_scaled))
          merged_columns[i].array = magerr_aper
        if merged_columns[i].name.startswith('%s_mag_aper_' % f):
          naper = merged_columns[i].name.split('_')[-1]
          naper = int(naper)
          flux_aper_n = c.__getitem__('%s_flux_aper_%d' % (f, naper))
          fluxerr_aper_n = c.__getitem__('%s_fluxerr_aper_%d' % (f, naper))
          fluxerr_aper_n_scaled = fluxerr_aper_n * skyrms_factor
          # sn_aper_n = flux_aper_n / fluxerr_aper_n
          sn_aper_n = flux_aper_n / fluxerr_aper_n_scaled
          mag_aper_n = np.where(sn_aper_n >= 1., c.__getitem__('%s_mag_aper_%d' % (f,naper)),
                                99.0)
          mag_aper_n = np.where(fluxerr_aper_n==0, -99.0, mag_aper_n)
          merged_columns[i].array = mag_aper_n
        if merged_columns[i].name.startswith('%s_magerr_aper_' % f):
          naper = merged_columns[i].name.split('_')[-1]
          naper = int(naper)
          flux_aper_n = c.__getitem__('%s_flux_aper_%d' % (f, naper))
          fluxerr_aper_n = c.__getitem__('%s_fluxerr_aper_%d' % (f, naper))
          fluxerr_aper_n_scaled = fluxerr_aper_n * skyrms_factor
          # sn_aper_n = flux_aper_n / fluxerr_aper_n
          sn_aper_n = flux_aper_n / fluxerr_aper_n_scaled
          # magerr_aper_n = np.where(fluxerr_aper_n==0, 0., 
          #                          c.__getitem__('%s_magerr_aper_%d' % (f,naper)))
          # magerr_aper_n = np.where(sn_aper_n>=1., magerr_aper_n, 
          #                          self.zeropoints[f]-2.5*np.log10(fluxerr_aper_n))
          magerr_aper_n = np.where(fluxerr_aper_n==0., 0., 
                                   2.5*np.log10(1.+1./sn_aper_n))
          magerr_aper_n = np.where(sn_aper>=1., magerr_aper, 
                                   self.zeropoints[f]-2.5*np.log10(fluxerr_aper_n_scaled))
          merged_columns[i].array = magerr_aper_n
      # OK, special treatment for mag_aper_1, which should be 0.4 arcsec 
      # diameter aperture
      for i in range(len(merged_columns)):
        if merged_columns[i].name == '%s_mag_aper_1' % f:
          j_mag_aper_1 = i
        elif merged_columns[i].name == '%s_magerr_aper_1' % f:
          j_magerr_aper_1 = i
        elif merged_columns[i].name == '%s_mag_aper_3' % f:
          j_mag_aper_3 = i
        elif merged_columns[i].name == '%s_magerr_aper_3' % f:
          j_magerr_aper_3 = i
      flux_aper_1 = c.__getitem__('%s_flux_aper_1' % f)
      fluxerr_aper_1 = c.__getitem__('%s_fluxerr_aper_1' % f)
      fluxerr_aper_1_scaled = fluxerr_aper_1 * skyrms_factor
      # sn_aper_1 = flux_aper_1 / fluxerr_aper_1
      sn_aper_1 = flux_aper_1 / fluxerr_aper_1_scaled
      mag_aper_1_1sig = np.where(sn_aper_1 >= 1., merged_columns[j_mag_aper_1].array,
                                 merged_columns[j_magerr_aper_1].array)
      merged_columns += [pyfits.Column(name='%s_mag_aper_1_1sig' % f,
                                       array=mag_aper_1_1sig, format='D')]
      flux_aper_3 = c.__getitem__('%s_flux_aper_3' % f)
      fluxerr_aper_3 = c.__getitem__('%s_fluxerr_aper_3' % f)
      fluxerr_aper_3_scaled = fluxerr_aper_3 * skyrms_factor
      # sn_aper_3 = flux_aper_3 / fluxerr_aper_3
      sn_aper_3 = flux_aper_3 / fluxerr_aper_3_scaled
      mag_aper_3_1sig = np.where(sn_aper_3 >= 1., merged_columns[j_mag_aper_3].array,
                                 merged_columns[j_magerr_aper_3].array)
      merged_columns += [pyfits.Column(name='%s_mag_aper_3_1sig' % f, 
                                       array=mag_aper_3_1sig, format='D')]
    print "len(merged_columns)", len(merged_columns)
    assert len(merged_columns) <= 999, "Warning: maximum number of columns allowed in a FITS table is 999!"
    coldefs = pyfits.ColDefs(merged_columns)
    tbhdu = pyfits.new_table(coldefs)
    merged_table = os.path.join(self.homedir, self.merged_table)
    print merged_table
    if os.path.exists(merged_table):
      os.remove(merged_table)
    tbhdu.writeto(merged_table)

class run_hst_photom(HSTPipeline):
  pass


class HSTPipeHot(HSTPipeline):
  """
  For hot-mode SExtractor run.
  """
  def __init__(self, cluster_name, paramfile_name):
    super(HSTPipeHot, self).__init__(cluster_name, paramfile_name)
    paramfile = os.path.join(surfsup_dir, paramfile_name)
    c = yaml.load(open(paramfile, 'rb'))
    # re-define output catalog name
    scale = c['scale'] * 1000  # in milli-arcsec
    scale = int(round(scale))
    scalestr = '%2dmas' % scale
    for b in self.filters:
      self.sex_catalog[b] = '%s_%s_%s_hot.cat' % (b, self.cluster_name.lower(), scalestr)
      self.sex_fitstable[b] = '%s_%s_%s_hot.fits' % (b, self.cluster_name.lower(), scalestr)
      if os.path.exists(self.sex_catalog[b]):
        os.remove(self.sex_catalog[b])
      if os.path.exists(self.sex_fitstable[b]):
        os.remove(self.sex_fitstable[b])
    # hot-mode-specific parameters
    self.detect_minarea = c['detect_minarea'][0]
    self.detect_thresh = c['detect_thresh'][0]
    self.clean_param = c['clean_param'][0]
    self.segimage = '%s_%s_%s_hot_seg.fits' % (self.detectband, self.cluster_name.lower(), scalestr)

  def run_sextractor(self, filters=[], sex_kwargs={}):
    if not sex_kwargs.has_key('detect_minarea'):
      sex_kwargs['detect_minarea'] = self.detect_minarea
    if not sex_kwargs.has_key('detect_thresh'):
      sex_kwargs['detect_thresh'] = self.detect_thresh
    if not sex_kwargs.has_key('clean_param'):
      sex_kwargs['clean_param'] = self.clean_param
    if not sex_kwargs.has_key('checkimage_name'):
      sex_kwargs['checkimage_name'] = self.segimage
    super(HSTPipeHot, self).run_sextractor(filters=filters, **sex_kwargs)

class HSTPipeCold(HSTPipeline):
  """
  For cold-mode SExtractor run.
  """
  def __init__(self, cluster_name, paramfile_name):
    super(HSTPipeCold, self).__init__(cluster_name, paramfile_name)
    paramfile = os.path.join(surfsup_dir, paramfile_name)
    c = yaml.load(open(paramfile, 'rb'))
    # re-define output catalog name
    scale = c['scale'] * 1000  # in milli-arcsec
    scale = int(round(scale))
    scalestr = '%2dmas' % scale
    for b in self.filters:
      self.sex_catalog[b] = '%s_%s_%s_cold.cat' % (b, self.cluster_name.lower(), scalestr)
      self.sex_fitstable[b] = '%s_%s_%s_cold.fits' % (b, self.cluster_name.lower(), scalestr)
      if os.path.exists(self.sex_catalog[b]):
        os.remove(self.sex_catalog[b])
      if os.path.exists(self.sex_fitstable[b]):
        os.remove(self.sex_fitstable[b])
    # cold-mode-specific parameters
    self.detect_minarea = c['detect_minarea'][1]
    self.detect_thresh = c['detect_thresh'][1]
    self.clean_param = c['clean_param'][1]
    self.segimage = '%s_%s_%s_cold_seg.fits' % (self.detectband, self.cluster_name.lower(), scalestr)

  def run_sextractor(self, filters=[], sex_kwargs={}):
    if not sex_kwargs.has_key('detect_minarea'):
      sex_kwargs['detect_minarea'] = self.detect_minarea
    if not sex_kwargs.has_key('detect_thresh'):
      sex_kwargs['detect_thresh'] = self.detect_thresh
    if not sex_kwargs.has_key('clean_param'):
      sex_kwargs['clean_param'] = self.clean_param
    if not sex_kwargs.has_key('checkimage_name'):
      sex_kwargs['checkimage_name'] = self.segimage
    super(HSTPipeCold, self).run_sextractor(filters=filters, **sex_kwargs)

class CombineColdHot(object):
  # combines the two modes
  def __init__(self, hotpipe, coldpipe, hotstart=10000):
    self.hotpipe = hotpipe
    self.coldpipe = coldpipe
    self.filters = self.hotpipe.filters
    self.objectid = np.zeros(0, 'int')
    self.hotstart = hotstart

  def hot_addition(self, kron_factor=3.0):
    # search in the hot-mode catalog and select objects that are at least 
    # kron_factor*kron_radius away from any cold-mode object
    # The objects added from hot mode will be their ID in hot mode + hotstart
    b0 = self.hotpipe.detectband
    chot = sextractor(self.hotpipe.sex_catalog[b0])
    ccold = sextractor(self.coldpipe.sex_catalog[b0])
    self.objectid = np.concatenate([self.objectid, ccold.number])
    hot_id = []
    for i in range(len(chot)):
      x = chot.x_image[i]
      y = chot.y_image[i]
      # calculate distances to all cold-mode objects
      dist = np.sqrt((x-ccold.x_image)**2 + (y-ccold.y_image)**2)
      min_cold_index = np.argsort(dist)[0]  
      # the index in cold catalog of the nearest neighbor
      mindist = dist.min()
      if mindist > (kron_factor * ccold.kron_radius[min_cold_index]):
        hot_id += [chot.number[i]]
    print "Number of hot additions: ", len(hot_id)
    new_hot_id = np.array(hot_id) + self.hotstart
    self.objectid = np.concatenate([self.objectid, new_hot_id])

  def combine_catalog(self):
    # combine the cold and hot mode catalogs.
    for b in self.filters:
      chot = sextractor(self.hotpipe.sex_catalog[b])
      ccold = sextractor(self.coldpipe.sex_catalog[b])
      output_name = self.hotpipe.sex_catalog[b].replace('hot','warm')
      f = open(output_name, 'wb')
      f.write(ccold._header)
      # include ALL objects from cold catalog
      for i in range(len(ccold)):
        f.write(ccold.line(i))
      # now incldue objects from hot catalog that are added by hot_addition
      for j in range(len(chot)):
        if (chot.number[j] + self.hotstart) in self.objectid:
          f.write(chot.line(j))
      f.close()
      print "Hot-cold catalog in %s complete." % b.upper()

  def combine_seg(self):
    # combine the hot and cold mode segmentation images
    hotseg = pyfits.getdata(self.hotpipe.segimage)
    coldseg = pyfits.getdata(self.coldpipe.segimage)
    # First, zero-out hot-mode objects not in final catalog
    b0 = self.hotpipe.detectband
    chot = sextractor(self.hotpipe.sex_catalog[b0])
    print "Zeroing out hot garbage..."
    for i in range(len(chot)):
      if (chot.number[i] + self.hotstart) not in self.objectid:
        hotseg = np.where(hotseg-chot.number[i]==0, 0, hotseg)
    # now combine the two segmentation images; replace coldseg with hotseg 
    # pixels only if coldseg==0
    print "Merging warm segmentation map..."
    warmseg = np.where((coldseg==0)&(hotseg>0), hotseg+self.hotstart, coldseg)
    # now write to an image 
    warmsegname = self.hotpipe.segimage.replace('hot', 'warm')
    if os.path.exists(warmsegname):
      os.remove(warmsegname)
    pyfits.append(warmsegname, warmseg, 
                  pyfits.getheader(self.coldpipe.segimage))

