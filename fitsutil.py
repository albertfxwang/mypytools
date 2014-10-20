#!/usr/bin/env python

from numpy import *
import pyfits   # requires pyfits 3.1 or later
from pygoods import *
import os

colformats = {'d':'J', 'f':'D', 's':'A'}
strformats = {'I':'%d','J':'%d','K':'%d','E':'%f','D':'%f','A':'%s'}

def sex2fits(c, fitsname, booleancols=[]):
	"""	usage: sex2fits(c, fitsname, booleancols=[])
	c -- input sextutils.sextractor catalog instance
	fitsname -- output file name of the binary FITS table
	booleancols -- column names that should be converted to boolean values"""
	fitscols = []
	# construct all the columns
	for i in range(len(c._d)):
		colname = c._colnames[i]
		coltype = c._type[colname]
		colfmt = colformats[coltype]
		if coltype == 's':
			slen = c._fmt[colname][1:-1]  # length of string
			colfmt = slen + colfmt
		colarray = c.__getattribute__(colname)
		# catch the Boolean array of 0 or 1
		if colname in booleancols:
			colfmt = 'L'
			colarray = where(colarray==0, False, True)  # convert to boolean array
		col = pyfits.Column(name=colname, format=colfmt, array=colarray)
		fitscols += [col]
	# create table header unit
	cols = pyfits.ColDefs(fitscols)
	tbhdu = pyfits.new_table(cols)
	hdu = pyfits.PrimaryHDU(array([]))  # create a primary HDU with an empty list
	thdulist = pyfits.HDUList([hdu, tbhdu])
	thdulist.writeto(fitsname)
	return 0
	
def fits2sex(fitsfile, sexfile, booleancols=[]):
	# writes the FITS table into the SExtractor catalog format
	h = pyfits.open(fitsfile)
	columns = h[1].columns
	colnames = []; formats = []
	nrow = shape(h[1].data)[0]  # number of rows (entries)
	ncol = len(columns)
	# read column information
	for i in range(ncol):
		colnames += [columns[i].name]
		formats += [columns[i].format]
	# write to sexfile
	f = open(sexfile, 'w')
	# first write headers
	for i in range(ncol):
		f.write('# %d %s\n' % ((i+1), colnames[i]))
	# now write contents
	for j in range(nrow):
		for i in range(ncol):
			if colnames[i] in booleancols:
				f.write('%d ' % (1 if h[1].data[j][i] else 0))
			f.write('%s ' % h[1].data[j][i])
		f.write('\n')
	f.close()
	return 0
	
	
def newfits_subset(oldfitsfile, newfitsfile, selection):
	# create a new fits table newfitsfile from a subset of oldfitsfile
	# selection is a boolean array telling us which entries from oldfitsfile to select
	c = Ftable(oldfitsfile)
	newcols = []
	for i in range(len(c._column_names)):
		colname = c._column_names[i]
		newarray = compress(selection, c.__getitem__(colname))
		col = pyfits.Column(name=colname, format=c.d.formats[i], array=newarray)
		newcols += [col]
	newcoldefs = pyfits.ColDefs(newcols)
	newhdu = pyfits.new_table(newcoldefs)
	newhdu.writeto(newfitsfile)
	return 0
	

def update_column(filename, colname, newarray):
	# update the column with a new array; the user should make sure that 
	# the length of the array is the same as before
	h = pyfits.open(filename, 'update')
	h[1].data.field(colname)[:] = newarray
	h.flush()
	h.close()

def replace_nan_num(filename, columns, value_dic):
	"""
	Replace the nan's with numbers. Give a list of column names and the
	corresponding values for replacement (as a dictionary). If a column name 
	is not found in value_dic, it is defaulted to zero.
	"""
	h = pyfits.open(filename, mode='update')
	for col in columns:
		if value_dic.has_key(col):
			val = value_dic[col]
		else:
			val = 0
		data = h[1].data.field(col)
		h[1].data.field(col)[:] = where(isnan(data), val, data)
	h.flush()
	h.close()


def change_column_names(filename, old_colnames, new_colnames):
	"""
	Change the name of a column.
	Pyfits does not really provide a convenient function to do this, so I'll have to 
	create a new table based on the old table, just use a different column name.
	"""
	os.system('mv %s %s.copy' % (filename, filename))
	c = pyfits.open(filename+".copy")
	tbhdu = c[1]
	ncol = len(tbhdu.data.columns)
	newcols = []
	for i in range(ncol):
		colname = tbhdu.data.columns[i].name
		colfmt = tbhdu.data.formats[i]
		colarr = tbhdu.data.field(colname)
		for j in range(len(old_colnames)):
			if tbhdu.data.columns[i].name == old_colnames[j]:
				colname = new_colnames[j]
				break
				#print colname
		newcols += [pyfits.Column(name=colname, format=colfmt, array=colarr)]
	newcols = pyfits.ColDefs(newcols)
	#print newcols
	newhdu = pyfits.new_table(newcols)
	newhdu.writeto(filename)
	c.close()
	os.system('rm %s.copy' % filename)
	
def add_columns(filename, colnames, newarrays, formats):
	# append new columns to an existing table
	oldfile = filename + ".OLD"
	try:
		os.system('mv %s %s' % (filename, oldfile))
	except:
		os.system('rm %s' % oldfile)
		os.system('mv %s %s' % (filename, oldfile))
	# write the new columns to a temporary table
	newcols = []
	for i in range(len(colnames)):
		newcols += [pyfits.Column(name=colnames[i],format=formats[i],
		                          array=newarrays[i])]
	newcols = pyfits.ColDefs(newcols)
	tbhdu = pyfits.new_table(newcols)
	try:
		tbhdu.writeto('temp.fits')
	except:
		os.system('rm temp.fits')
		tbhdu.writeto('temp.fits')
	# Now read in the old table and merge
	h1 = pyfits.open(oldfile)
	h2 = pyfits.open('temp.fits')
	h = h1[1].columns + h2[1].columns  # merge the columns
	newhdu = pyfits.new_table(h)
	newhdu.writeto(filename)
	return 0
	
def add_comment_from_file(fitsfile, commentfile, cmttag='#'):
	# add comments from a text file to the header of fitsfile
	# cmttag is the tag identifying comment lines in commentfile
	# all the comments are under the "COMMENT" keyword within the header
	f = open(commentfile)
	lines = f.readlines()
	commentlines = []
	# gather all the comment lines
	for l in lines:
		if l.startswith(cmttag):
			commentlines += [l.strip()]
	# now write to the new header
	h = pyfits.open(fitsfile, 'update')
	for cl in commentlines:
		h[1].header['COMMENT'] = cl
	h.flush()
	h.close()
	
def merge_tables(table1, table2, newtabname, mode='left'):
	"""	usage: merge_tables(table1, table2, newtabname, mode='left')
	table1, table2 are two fits tables (Ftable instances) to be merged.
	mode:
	'left': only merge columns that are in table1; ignore other columns in table2.
	'right': only merge columns that are in table2; ignore other columns in table1.
	"""
	columns = []
	if mode == 'left':
		for i in range(len(table1.Columns)):
			col = table1.Columns[i]
			if col in table2.Columns:  # if col is also in table2
				a1 = table1.__getitem__(col)
				a2 = table2.__getitem__(col)
				a_all = concatenate([a1,a2])
				fmt = table1.d.formats[i]
				columns += [pyfits.Column(name=col, format=fmt, array=a_all)]
	elif mode == 'right':
		for i in range(len(table2.Columns)):
			col = table2.Columns[i]
			if col in table1.Columns:  # if col is also in table1
				a2 = table2.__getitem__(col)
				a1 = table1.__getitem__(col)
				a_all = concatenate([a1,a2])
				fmt = table2.d.formats[i]
				columns += [pyfits.Column(name=col, format=fmt, array=a_all)]
	cols = pyfits.ColDefs(columns)
	tbhdu = pyfits.new_table(cols)
	tbhdu.writeto(newtabname)
		