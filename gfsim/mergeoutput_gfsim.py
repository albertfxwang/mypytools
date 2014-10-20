#!/usr/bin/env python

import numpy as N
from pygoods import Ftable, sextractor
import fitsutil
import pyfits
import os, glob
import time

"""
Merge all the gfit*_matched.fits; this task simply merges all the FITS tables
and does not do any kind of matching at all.
"""

def mergeoutput_gfsim(root,bands,output_root='merged_gfsim',
                      exclude_columns=['filename']):
	""" 
	Make one merged catalog for each band.
	exclude_columns: a kludge to skip certain columns that do not exist in all
	runs. This should not have happened in the first place...
	"""
	timestr = time.strftime('%y%m%d')
	for b in bands:
		broot = root+'_'+b   # the directory of the given band
		mcats = glob.glob(broot+'/gfit*_matched.fits')
		# print "mcats[0]", mcats[0]
		c0 = Ftable(mcats[0])
		colnames = c0.Columns
		# colnames = c0._colnames
		for col in colnames:
			if col in exclude_columns:
				colnames.remove(col)
		colformats = c0.d.formats
		colarrays_dic = {}
		# initialize colarrays
		for i in range(len(colnames)):
			col = colnames[i]
			colarrays_dic[col] = getattr(c0,col)
		# Now concatenate columns from all catalogs
		for i in range(1,len(mcats)):   # we already added the first one
			ci = Ftable(mcats[i])
			# print mcats[i]
			for col in colnames:
				# print len(colarrays_dic[col]), len(getattr(ci, col)), col
				colarrays_dic[col] = N.concatenate([colarrays_dic[col],getattr(ci,col)])
		# Now write the output catalog
		outname = output_root+'_'+b+'_'+timestr+'.fits'
		columns = []
		for i in range(len(colnames)):
			columns += [pyfits.Column(name=colnames[i],format=colformats[i],
			                          array=colarrays_dic[colnames[i]])]
		ColDefs = pyfits.ColDefs(columns)
		tbhdu = pyfits.new_table(ColDefs)
		tbhdu.writeto(outname)
		print "%s done." % outname


def update_sex_colnames(catalog, oldnames, newnames):
	# Update the SExtractor catalog headers
	assert len(oldnames) == len(newnames)
	f = open(catalog)
	lines = f.readlines()
	f.close()
	# loop through all lines and look for match from oldnames
	for i in range(len(lines)):
		if lines[i][0] == '#':
			for j in range(len(oldnames)):
				if oldnames[j] in lines[i]:
					lines[i] = lines[i].replace(oldnames[j], newnames[j])
	f2 = open(catalog, 'wb')
	for i in range(len(lines)):
		f2.write(lines[i])
	f2.close()

