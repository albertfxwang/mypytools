#!/usr/bin/env python
##==========================================================================##
## For runs using sexsim/sextractor_sim.py, use the collect_results method.
##==========================================================================##
 

from numpy import *
import os,glob
from pygoods import *
import pyfits
#import simutil as su


zeropoints = {}
zeropoints['F435W'] = 25.65288
zeropoints['F606W'] = 26.49341
zeropoints['F775W'] = 25.64053
zeropoints['F850LP'] = 24.84315
zeropoints['F098M'] = 26.270
zeropoints['F105W'] = 26.261
zeropoints['F125W'] = 26.237
zeropoints['F160W'] = 25.932
zeropoints['JH'] = 26.057

def_inputcols = {'x_in':1,'y_in':2,'mag_in':3,'gtype':4,'Re_in':5,'Ell_in':6,'theta_in':7,
	'z_in':14,'ebmv_in':15,'M1500_in':16,'U_mag_in':9,'i_mag_in':10,'z_mag_in':11,'Y_mag_in':12,
	'J_mag_in':13}

def_inputcols2 = {'x_in':1,'y_in':2,'mag_in':3,'gtype':4,'Re_in':5,'Ell_in':6,'theta_in':7}

def_inputcols_templates = {'x_in':1, 'y_in':2, 'mag_in':3, 'template':4, 'tempsize':5, 
	'ell_in':6, 'pa_in':7, 'magerr_in':9, 'xtemp':10, 'ytemp':11, 'r_h_in':12,
	'multi':13}

def flux2mag(flux,zpt):
   Fneg = 10.**((zpt - 99.0) / 2.5)
   try:
      flux = float(flux)
   except:
      flux = maximum(flux, Fneg)
      return zpt - 2.5 * log10(flux)
   else:
      if flux <= 0:
         return 99.0
      else:
         return zpt - 2.5 * log10(flux)
   
def mergeoutput(dirs, root='run*m', bands=['H','J'],
   inputband='H',inputcols_dic=def_inputcols,inputcols_int=['gtype'],paramfile=None):
   """	usage: mergeoutput(dirs, root='run*m', bands=['H','B','V','I','Z','Y','J'],
   inputband='H',inputcols_dic=def_inputcols,inputcols_int=['gtyep'], naper=1)
   dirs -- a list of all simulation directories
   root -- the file name of the *.sim file without the extension
   bands -- all the bands in the simulations
   inputband -- the band from which the input parameters are taken from
   inputcols_dic -- a dictionary of (key, value) pairs indicating the column names (key)
   	of the nth column (value) in the input catalog 
   inputcols_int -- the column names whose values should be integer and not float
   #naper -- number of apertures used in the SExtractor run   
   """
   print inputcols_dic
   runs = []
   for d in dirs:
      if d[-1] == '/': d = d[:-1]
      runs += glob.glob(d+'/'+root+'_%s' % inputband)
   print 'len(runs)', len(runs)
   cats = {}  # a dictionary of SExtractor catalogs
   gallist = {}  # input catalog
   BIG_DICT = {}   
   BIG_DICT['detect'] = zeros(0, 'bool')
   BIG_DICT['mag_in'] = {}
   ALLCOLUMNS = ['detect']
   ALLFORMATS = ['L']
   
   
   # the following input parameters taken from detection-band catalog
   inputcols = {}
   for k in inputcols_dic.keys():
      if k in inputcols_int:
         #inputcols[k] = zeros(0, 'int')
         BIG_DICT[k] = zeros(0, 'int')
      elif k != 'mag_in': # initialize input magnitudes separately
         #inputcols[k] = zeros(0)
         BIG_DICT[k] = zeros(0, 'float')
      elif k == 'mag_in':
      	for b in bands:
      		BIG_DICT['mag_in'][b] = zeros(0, 'float')
      if k not in ALLCOLUMNS:
         ALLCOLUMNS += [k]
         if k in inputcols_int:
            ALLFORMATS += ['I']
         else:
            ALLFORMATS += ['D']
   
   firstcat = 1
   for run in runs:
      subruns = glob.glob(run+'/run*.cat')
      for srun in subruns:
         i = srun.split('/')[-1]
         i = i.split('.')[0]
         i = i[3:]; i = int(i)
         f = open(srun)
         lines = f.readlines()
         f.close()
         ndetobj = 0
         for l in lines:
            if l[0] != '#':
               ndetobj += 1
         #if len(lines) <= 59:
         if ndetobj == 0:
            print "no detected objects; something wrong with this run"
            continue
         while firstcat:
            if paramfile == None:
               # determine the columns if paramfile == None
               c = sextractor(srun)
               ncol = 0
               OUTPUT_COLUMNS = c._colnames
               OUTPUT_FORMATS = c._type
               for k in OUTPUT_FORMATS.keys():
                  if k not in ALLCOLUMNS:
                     if 'vector_assoc' not in k:
                        ALLCOLUMNS += [k]
                        if OUTPUT_FORMATS[k] == 'd':
                           ALLFORMATS += ['I']
                        else:
                           ALLFORMATS += ['D']
                  if OUTPUT_FORMATS[k] == 'd':
                     OUTPUT_FORMATS[k] = 'int'
                  else:
                     OUTPUT_FORMATS[k] = 'float'            
               for colname in OUTPUT_COLUMNS:
                  if 'vector_assoc' not in colname: 
                     BIG_DICT[colname] = {}
                     for b in bands:
                        # avoid vector_assocs; match with input columns instead
                        BIG_DICT[colname][b] = zeros(0, OUTPUT_FORMATS[colname])
                     ncol += 1
               firstcat = 0
               print "%d columns from output catalogs included." % ncol
            else:
               # read the columns from a file
               # each line of the file should contain column name and data format 
               # (int, float, or bool)
               OUTPUT_COLUMNS = []
               OUTPUT_FORMATS = {}
               ncol = 0
               f = open(paramfile)
               lines = f.readlines()
               for line in lines:
               	if (line[0] != '#') or (line.strip() != ""):
                     l = line.split()
                     colname = l[0].lower()
                     colformat = l[1]  # 'int', 'float', or 'bool'
                     if colname[-1] == ')':
                        # there are more than one columns associated with this name
                        s = colname[:-1].split('(')
                        numcol = int(s[1])
                        colname_base = s[0]
                        if colname_base not in ALLCOLUMNS:
                           for j in range(numcol):
                              if j == 0:
                                 ALLCOLUMNS += [colname_base]
                                 OUTPUT_COLUMNS += [colname_base]
                                 OUTPUT_FORMATS[colname_base] = colformat
                              else:
                                 ALLCOLUMNS += [colname_base+'_%d'%j]
                                 OUTPUT_COLUMNS += [colname_base+'_%d'%j]
                                 OUTPUT_FORMATS[colname_base+'_%d'%j] = colformat
                              if colformat == 'int':
                                 ALLFORMATS += ['I']
                              elif colformat == 'float':
                                 ALLFORMATS += ['D']
                              else:
                                 ALLFORMATS += ['L']
                     
                     elif colname not in ALLCOLUMNS:
                        ALLCOLUMNS += [colname]
                        OUTPUT_COLUMNS += [colname]
                        OUTPUT_FORMATS[colname] = colformat
                        if colformat == 'int':
                           ALLFORMATS += ['I']
                        elif colformat == 'float':
                           ALLFORMATS += ['D']
                        else:
                           ALLFORMATS += ['L']
               for colname in OUTPUT_COLUMNS:
                  if 'vector_assoc' not in colname: 
                     BIG_DICT[colname] = {}
                     for b in bands:
                        # avoid vector_assocs; match with input columns instead
                        BIG_DICT[colname][b] = zeros(0, OUTPUT_FORMATS[colname])
                     ncol += 1
               firstcat = 0
               print "%d columns from output catalogs included." % ncol
         #runroot = srun.split('/')[0]
         runroot = run[:-2]  # e.g. run1m
         for b in bands:
            cats[b] = sextractor(runroot+'_%s/run%d.cat' % (b, i))
            # the SExtractor output file
            if len(b) == 1:  # if detection is using a single band
               # gallist[b] = sextractor(runroot+'_%s/gal%d.list' % (b, i))
               gallist[b] = sextractor(runroot+'_%s/glart%d.list' % (b, i))
               # the input file for the given run
            else: 
               # a composite detection band does not have input catalog, since the image is 
               # calculated from two (or more) bands
               gallist[b] = None
         ninput = len(gallist[inputband])  # number of input galaxies for each run (e.g. 40)
         detect_srun = zeros(ninput, 'bool')
         for b in bands:
            BIG_DICT_BAND = {}
            # record input magnitudes
            BIG_DICT_BAND['mag_in'] = {}
            BIG_DICT_BAND['mag_in'][b] = gallist[b].__getattribute__('_%d'%inputcols_dic['mag_in'])
            # First, get the input values
            if b == inputband:  # band-independent attributes
               # first collect the input parameters and attach back to BIG_DICT
               for kin in inputcols_dic.keys():
                  if kin == 'gtype':
                     gtype0 = array(gallist[b].__getattribute__('_%d'%inputcols_dic['gtype']))
                     gtype0 = where(gtype0 == 'expdisk', 0, 1).astype('int')  # convert ('expdisk', 'devauc') to (0, 1)
                     BIG_DICT['gtype'] = concatenate((BIG_DICT['gtype'], gtype0))
                  elif kin != 'mag_in':
                     BIG_DICT[kin] = concatenate((BIG_DICT[kin], gallist[b].__getattribute__('_%d'%inputcols_dic[kin])))
               #detect_b = zeros(ninput, 'int')
            # initialize the arrays for each band
       
            for colname in OUTPUT_COLUMNS:
               BIG_DICT_BAND[colname] = ones(ninput, OUTPUT_FORMATS[colname]) * -1
            
            # now read in the detected sources
            # Now match the input objects to the output catalog, if they are detected
            for j in range(ninput):
            	# if the object is detected...
               if (gallist[inputband]._1[j] in cats[b].vector_assoc) &\
                  (gallist[inputband]._2[j] in cats[b].vector_assoc_1):
                  ix = where(cats[b].vector_assoc==gallist[inputband]._1[j])[0]
                  iy = where(cats[b].vector_assoc_1==gallist[inputband]._2[j])[0]
                  ii = intersect1d(ix,iy)[0]
                  if b == inputband: detect_srun[j] = True
                  for colname in OUTPUT_COLUMNS:
                  	# paste the values read from the output catalog into BIG_DICT_BAND
                  	BIG_DICT_BAND[colname][j] = cats[b].__getattribute__(colname)[ii]
                  
            # attach the new run at the end of each attribute
            for colname in OUTPUT_COLUMNS:
               if 'vector_assoc' not in colname:
            	   BIG_DICT[colname][b] = concatenate((BIG_DICT[colname][b], BIG_DICT_BAND[colname]))
            BIG_DICT['mag_in'][b] = concatenate((BIG_DICT['mag_in'][b], 
               BIG_DICT_BAND['mag_in'][b]))
            
         BIG_DICT['detect'] = concatenate((BIG_DICT['detect'], detect_srun))
   
   return BIG_DICT, ALLCOLUMNS, ALLFORMATS

def mergecat_txt(outcat, dirs, root='run*m_udf_z8', bands=['H','B','V','I','Z','Y','J'],\
   inputband='H', inputcols_dic=def_inputcols, inputcols_int=[], inputcols_str=[],
   paramfile=None):
   md, allcolumns, allformats = mergeoutput(dirs, root=root, bands=bands, inputband=inputband,
      inputcols_dic=inputcols_dic, inputcols_int=inputcols_int,paramfile=paramfile)
   # First build the headers
   header = ""
   header += "# 1 NUMBER\n"
   nc = 2
   for i in range(len(allcolumns)):
      colname = allcolumns[i]
      format = allformats[colname]
      if type(md[colname]) == type({}):
         for b in bands:
            header += "# %d %s_%s\n" % (nc, b, colname)
            nc += 1
      else:
         header += "# %d %s\n" % (nc, colname)
         nc += 1   
   
   
   # write the catalog
   f = open(outcat,'w')
   f.write(header)
   for i in range(len(md['detect'])):
      # iterate through all objects
      f.write('%d ' % i)  # write NUMBER
      for j in range(len(allcolumns)):
         colname = allcolumns[i]
         if type(md[colname]) == type({}):  # this column has different values for each band
            for b in bands:
               # if this column has integer values
               if allformats[colname] in ['I','J','L']:
                  f.write('%d ' % md[colname][b])
               elif allformats[colname] in ['D','E']:
                  f.write('%f ' % md[colname][b])
               else:
                  f.write('%s ' % md[colname][b])
         else:
            # this column has the same value for all bands
            if allformats[colname] in ['I','J','L']:
               f.write('%d ' % md[colname])
            elif allformats[colname] in ['D','E']:
               f.write('%f ' % md[colname])
            else:
               f.write('%s ' % md[colname])
      f.write('\n')

   f.close()

def mergecat_fits(outcat, dirs, root='run*m_udf_z8', bands=['H','B','V','I','Z','Y','J'],\
   inputband='H',inputcols_dic=def_inputcols,inputcols_int=['gtype'],paramfile=None):
   # collects the results from all the simulation catalogs and write to a *FITS* table
   # collect the results
   if len(glob.glob(outcat)) > 0:
   	raise ValueError, "%s already exists." % outcat
   md, allcolumns, allformats = mergeoutput(dirs, root=root, bands=bands, inputband=inputband, inputcols_dic=inputcols_dic, 
      inputcols_int=inputcols_int,paramfile=paramfile)
   # now build the FITS table
   fitscols = []
   for i in range(len(allcolumns)):
      colname = allcolumns[i]
      format = allformats[i]
      if type(md[colname]) == type(array([])):
         fitscols += [pyfits.Column(name=colname, format=format, array=md[colname])]
      else:
         for b in bands:
            bcolname = b + '_' + colname
            fitscols += [pyfits.Column(name=bcolname, format=format, 
               array=md[colname][b])]
   fitscols += [pyfits.Column(name='number', format='K', 
      array=arange(0, len(md['detect'])))]
   cols = pyfits.ColDefs(fitscols)
   tbhdu = pyfits.new_table(cols)
   hdu = pyfits.PrimaryHDU(array([]))
   thdulist = pyfits.HDUList([hdu, tbhdu])
   thdulist.writeto(outcat)
   	
   
def add_header(header, lastnum, newattr, newcomments=""):
   nh = lastnum + 1
   header = header + "#  %d %s    %s\n" % (nh, newattr, newcomments)
   return header, nh

def merge_skysim(outcat, dirs,root, bands=['JH','B','V','I','Z','Y','J','H']):
   # mergeout the outputs of SExtractor simulation on pure sky-noise image
   headers = ""
   headers += """
#   1 NUMBER          The number of run (e.g. 0 for run0.cat)
#   2 X_IMAGE         Object position along x                         [pixel]
#   3 Y_IMAGE         Object position along y                         [pixel]
#   4 ALPHA_J2000     Right ascension of barycenter (J2000)           [deg]
#   5 DELTA_J2000     Declination of barycenter (J2000)               [deg]
#   6 ISOAREAF_IMAGE  Isophotal area (filtered) above Detection thres [pixel**2]
#   7 ISOAREA_IMAGE   Isophotal area above Analysis threshold         [pixel**2]
#   8 THETA_IMAGE     Position angle (CCW/x)                          [deg]
#   9 ELLIPTICITY     1 - B_IMAGE/A_IMAGE
#  10 CLASS_STAR      S/G classifier output
"""
   nh = 10
   for b in bands:
      headers, nh = add_header(headers, nh, '%s_KRON_RADIUS' % b)
      headers, nh = add_header(headers, nh, '%s_FLUX_RADIUS' % b)
      nh = nh + 2 
      headers, nh = add_header(headers, nh, '%s_FWHM_IMAGE' % b)
      headers, nh = add_header(headers, nh, '%s_FLAGS' % b)
      headers, nh = add_header(headers, nh, '%s_IMAFLAGS_ISO' % b)
      headers, nh = add_header(headers, nh, '%s_NIMAFLAGS_ISO' % b)
      headers, nh = add_header(headers, nh, '%s_BACKGROUND' % b)
      headers, nh = add_header(headers, nh, '%s_MAG_ISO' % b)
      headers, nh = add_header(headers, nh, '%s_MAGERR_ISO' % b)
      headers, nh = add_header(headers, nh, '%s_FLUX_ISO' % b)
      headers, nh = add_header(headers, nh, '%s_FLUXERR_ISO' % b)
      headers, nh = add_header(headers, nh, '%s_MAG_AUTO' % b)
      headers, nh = add_header(headers, nh, '%s_MAGERR_AUTO' % b)
      headers, nh = add_header(headers, nh, '%s_FLUX_AUTO' % b)
      headers, nh = add_header(headers, nh, '%s_FLUXERR_AUTO' % b)
   broot_det = root+'_%s' % detectband
   runcats = glob.glob('%s/run*.cat' % broot_det)
   f = open(outcat, 'w')
   f.write(headers)
   for rc in runcats:
      sexcats = {}
      try:
         c0 = sextractor(rc)
      except:
         print "can't read %s" % rc
         continue
      cat = rc.split('/')[-1]
      rc = cat[3:-4]; rc = int(rc)  
      # write detection-band only attributes
      sexcats = {}
      for b in bands:
         broot = root + '_%s' % b
         sexcats[b] = sextractor(broot+'/run%d.cat' % rc)
      for i in range(len(c0)):
         f.write('%d %f %f ' % (rc, c0.x_image[i], c0.y_image[i]))
         f.write('%f %f ' % (c0.alpha_j2000[i], c0.delta_j2000[i]))
         f.write('%f %f ' % (c0.isoareaf_image[i], c0.isoarea_image[i]))
         f.write('%f %f %f ' % (c0.theta_image[i], c0.ellipticity[i], c0.class_star[i]))
         
         #broot = root + '_%s' % b
         #cb = sextractor(broot+'/run%d.cat' % rc)
         for b in bands:
            
            cb = sexcats[b]
            f.write('%f %f %f %f ' % (cb.kron_radius[i], cb.flux_radius[i], cb.flux_radius_1[i],
               cb.flux_radius_2[i]))
            f.write('%f %d %d %d ' % (cb.fwhm_image[i], cb.flags[i], cb.imaflags_iso[i],
               cb.nimaflags_iso[i]))
            f.write('%f %f %f ' % (cb.background[i], cb.mag_iso[i], cb.magerr_iso[i]))
            f.write('%f %f ' % (cb.flux_iso[i], cb.fluxerr_iso[i]))
            f.write('%f %f ' % (cb.mag_auto[i], cb.magerr_auto[i]))
            f.write('%f %f ' % (cb.flux_auto[i], cb.fluxerr_auto[i]))
         f.write('\n')
   f.close()

