#!/usr/bin/env python

from numpy import *
import os,glob
from pygoods import *
import pyfits
import string
#import simutil as su


zeropoints = {}
zeropoints['B'] = 25.65288
zeropoints['V'] = 26.49341
zeropoints['I'] = 25.64053
zeropoints['Z'] = 24.84315
zeropoints['Y'] = 26.261
zeropoints['J'] = 26.237
zeropoints['H'] = 25.932
zeropoints['JH'] = 26.057

def_inputcols = {'x_in':1,'y_in':2,'mag_in':3,'gtype':4,'Re_in':5,'Ell_in':6,'theta_in':7}

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

def getcolumns(dir, root, bands):
	"""
	Get the columns in all bands.
	"""
	OUTPUT_COLUMNS = []
	OUTPUT_FORMATS = {}
	ALLFORMATS = []
	for b in bands:
		runs = glob.glob(dir+'/'+root+'_'+b+'/run*.cat')
		print runs[0]
		c = sextractor(runs[0])
		for cn in c._colnames:
			bcn = b+'_'+cn  # prepend the band name to each column
			if ('vector_assoc' not in bcn) & ('number' not in bcn):
				# the returned list of column names will NOT include all the VECTOR_ASSOC_*
				# and NUMBER columns
				OUTPUT_COLUMNS += [bcn]
				OUTPUT_FORMATS[bcn] = c._type[cn]
				if c._type[cn] == 'd':
					ALLFORMATS += ['I']
				else:
					ALLFORMATS += ['D']
				if OUTPUT_FORMATS[bcn] == 'd':
					OUTPUT_FORMATS[bcn] = 'int'
				else:
					OUTPUT_FORMATS[bcn] = 'float'
	OUTPUT_COLUMNS += ['idsim']
	OUTPUT_FORMATS['idsim'] = 'int'
	ALLFORMATS += ['L']
	print "%d columns from output catalogs included." % len(OUTPUT_COLUMNS)
	return OUTPUT_COLUMNS, OUTPUT_FORMATS, ALLFORMATS
	
   
def mergeoutput(dir, root='multires_candels', bands=['f160w','ch1'],
   inputband='f160w',inputcols_dic=def_inputcols,
   nstart_id_sim=100000,ngal_iter=200,hiresbands=['f160w'],idcolnames=['number','objectid']):
   """	usage: mergeoutput(dirs, root='run*m', bands=['H','B','V','I','Z','Y','J'],
   inputband='H',inputcols_dic=def_inputcols, naper=1)
   dir -- the simulation directory
   root -- the file name of the *.sim file without the extension
   bands -- all the bands in the simulations
   inputband -- the band from which the input parameters are taken from
   inputcols_dic -- a dictionary of (key, value) pairs indicating the column names (key)
   	of the nth column (value) in the input catalog 
   #naper -- number of apertures used in the SExtractor run   
   NOTE: no need to include objects not detected in the detection band?
   """
   print inputcols_dic
   colname_dic = {}
   for i in range(len(bands)):
      colname_dic[bands[i]] = idcolnames[i]
   ALLCOLUMNS, OUTPUT_FORMATS, ALLFORMATS = getcolumns(dir, root, bands)
   # get all the columns from output catalogs
   BIG_DICT = {}
   BIG_DICT['detect'] = zeros(0, 'bool')
   BIG_DICT['mag_in'] = {}
   ALLCOLUMNS += ['detect']
   ALLFORMATS += ['L']
   rootdir = dir+'/'+root
   
   # the following input parameters taken from detection-band catalog
   inputcols = {}
   for k in inputcols_dic.keys():
      if k == 'gtype':
         #inputcols[k] = zeros(0, 'int')
         BIG_DICT[k] = zeros(0, 'int')
      elif k != 'mag_in': # initialize input magnitudes separately
         #inputcols[k] = zeros(0)
         BIG_DICT[k] = zeros(0, 'float')
      else:  # take input magnitude from each band, not just the input band
      	for b in bands:
      		BIG_DICT['mag_in'][b] = zeros(0, 'float')
      if k not in ALLCOLUMNS:
         ALLCOLUMNS += [k]
         if k == 'gtype':
            ALLFORMATS += ['I']
         else:
            ALLFORMATS += ['D']
   
   # create key-value pairs in BIG_DICT for all columns
   for k in ALLCOLUMNS:
   	if k not in BIG_DICT.keys():
   		BIG_DICT[k] = zeros(0, OUTPUT_FORMATS[k])
   		# will include the column "detect"
            
   # Now start matching the entries in the output catalogs...
   # The IDs in TFIT catalogs are calculated from the IDs in SExtractor catalogs by
   # id_tfit = nstart_id_sim + ngal_iter * niter + id_sex
   runs = glob.glob(rootdir+'_'+inputband+'/run*.cat')  # collect all the runs in inputband
   print "len(runs)", len(runs)
   for r in runs:
      niter = r.split('/')[-1]
      niter = niter.split('.')[0]
      niter = niter.split('_')[0]
      niter = niter[3:]; niter = int(niter)
      catexists = []
      for b in bands:
         catexists += [os.path.exists(rootdir+'_'+b+'/run%d_%s.cat'%(niter,b))]
      catexists = array(catexists)
      if not catexists.all():
         continue
      id_offset = nstart_id_sim + niter * ngal_iter
      BIG_DICT_RUN = {}
      BIG_DICT_RUN['mag_in'] = {}
      c_in = sextractor(rootdir+'_'+inputband+'/run%d_%s.cat'%(niter,inputband))
      gl_in = sextractor(rootdir+'_'+inputband+'/gal%d_%s.list'%(niter,inputband))
      BIG_DICT_RUN['idsim'] = ones(len(gl_in), 'int') * -1  # initialize the array for ID
      BIG_DICT_RUN['detect'] = zeros(len(gl_in), 'int')
      # the SExtractor output catalog in the input band... will use later
      for b in bands:
      	brootdir = rootdir+'_'+b
      	idcolname = colname_dic[b]
      	gl = sextractor(brootdir+'/gal%d_%s.list' % (niter,b))  # the input list 
      	if b == inputband:  # fetch the input values except for magnitude
      		for k in inputcols_dic.keys():
      			if k == 'gtype':
      				karr = gl.__getattribute__('_%d'%inputcols_dic[k])
      				BIG_DICT_RUN[k] = where(karr=='devauc', 1, 0)
      			elif k != 'mag_in':
      				BIG_DICT_RUN[k] = gl.__getattribute__('_%d'%inputcols_dic[k])
      	# now get the input magnitude for this band
      	BIG_DICT_RUN['mag_in'][b] = gl.__getattribute__('_%d'%inputcols_dic['mag_in']) 
      	brootdir = rootdir + '_' + b
      	c = sextractor(brootdir+'/run%d_%s.cat' % (niter, b))
      	# Now find the matches of input sources in the output catalog 
      	# first initialize all the columns in BIG_DICT_RUN for this band
      	for col in ALLCOLUMNS:
      		if col.startswith(b) & ('mag_in' not in col):
      			# don't initialize input magnitude... it is already done
      			BIG_DICT_RUN[col] = ones(len(gl), OUTPUT_FORMATS[col]) * -1
      	# Now figure out if this band is run with SExtractor or TFIT
      	if b in hiresbands:
      		# in hi-res band... match by input position
      		# this will be slow because it's looping through all the input objects...
      		# but I can't think of a more clever way of doing it so far
      		for i in range(len(gl)):
      			x = (c.vector_assoc==gl._1[i]) & (c.vector_assoc_1==gl._2[i])
      			if sum(x) > 0:  # if found a match
      				#print "match"
      				idx = arange(len(x))[x==True]
      				# calculate the simulation ID for this found object
      				num = c.number[idx]
      				idsim = nstart_id_sim + (niter-1) * ngal_iter + num
      				BIG_DICT_RUN['idsim'][i] = idsim
      				BIG_DICT_RUN['detect'][i] = 1
      				for col in ALLCOLUMNS:
      					if col.startswith(b) & ('mag_in' not in col) & ('detect' not in col):
      						scol = col[col.find('_')+1:]  # strip the column name of band name
      						BIG_DICT_RUN[col][i] = c.__getattribute__(scol)[idx]
      						# copy the value of the column scol found for the object
      	else:
      		# this is a low-res band run by TFIT... match by ID number
      		# match to the idsim taken from a high-res catalog
      		for i in range(len(gl)):
      			idsim = BIG_DICT_RUN['idsim'][i]
      			idcolumn = c.__getattribute__(idcolname)
      			if idsim in idcolumn:
      				idy = arange(len(c))[idcolumn==idsim][0]
      				for col in ALLCOLUMNS:
      					if col.startswith(b) & ('mag_in' not in col) & ('detect' not in col):
      						scol = col[col.find('_')+1:]  # strip the column name of band name
      						try:
      							BIG_DICT_RUN[col][i] = c.__getattribute__(scol)[idy]
      							# copy the value of the column scol found for the object
      						except:
      							# if for some reason col is not one of the columns in c, for this run
      							BIG_DICT_RUN[col][i] = -1
      # OK... now with BIG_DICT_RUN, paste those values back to BIG_DICT
      for k in BIG_DICT_RUN.keys():
      	if k == 'mag_in':
      		for b in bands:
      			BIG_DICT[k][b] = concatenate((BIG_DICT[k][b], BIG_DICT_RUN[k][b]))
      	else:
      		BIG_DICT[k] = concatenate((BIG_DICT[k], BIG_DICT_RUN[k]))
   
   return BIG_DICT, ALLCOLUMNS, ALLFORMATS

def mergecat_txt(outcat, dirs, root='run*m_udf_z8', bands=['H','B','V','I','Z','Y','J'],\
   inputband='H', inputcols_dic=def_inputcols):
   md, allcolumns, allformats = mergeoutput(dirs, root=root, bands=bands, inputband=inputband,
      inputcols_dic=inputcols_dic)
   # First build the headers
   header = ""
   nc = 1
   for i in range(len(allcolumns)):
      colname = allcolumns[i]
      format = allformats[i]
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
      #f.write('%d ' % i)  # write NUMBER
      for j in range(len(allcolumns)):
         colname = allcolumns[j]
         if type(md[colname]) == type({}):  # this column has different values for each band
            for b in bands:
               # if this column has integer values
               if allformats[j] in ['I','J','L']:
                  f.write('%d ' % md[colname][b][i])
               elif allformats[j] in ['D','E']:
                  f.write('%f ' % md[colname][b][i])
               else:
                  f.write('%s ' % md[colname][b][i])
         else:
            # this column has the same value for all bands
            if allformats[j] in ['I','J','L']:
               f.write('%d ' % md[colname][i])
            elif allformats[j] in ['D','E']:
               f.write('%f ' % md[colname][i])
            else:
               f.write('%s ' % md[colname][i])
      f.write('\n')

   f.close()

def mergecat_fits(outcat, dirs, root='uvimos_goodss', bands=['f435w','uvimos'],
   inputband='f435w', inputcols_dic=def_inputcols,
   idcolnames=['number','objectid'], hiresbands=['f435w']):
   # collects the results from all the simulation catalogs and write to a *FITS* table
   # collect the results
   if len(glob.glob(outcat)) > 0:
   	raise ValueError, "%s already exists." % outcat
   if len(idcolnames)==0:  
   	idcolnames = ['number'] * len(bands)  # the default column name for the ID column in the catalogs
   md, allcolumns, allformats = mergeoutput(dirs, root=root, bands=bands, 
                                inputband=inputband, 
                                inputcols_dic=inputcols_dic,
                                idcolnames=idcolnames, 
                                hiresbands=hiresbands)
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

