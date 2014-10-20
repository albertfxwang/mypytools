#!/usr/bin/env python

import numpy as N
import pyfits
from pygoods import sextractor, Ftable
import os, glob, cPickle, copy, subprocess
from pyraf import iraf
try:
	import pywcs
	WCS = pywcs.WCS
except:
	import astropy
	WCS = astropy.wcs.WCS

c = Ftable('gs_all_tf_h_130511a_multi_photoz_v1.0.fits')
rootdir = '/astro/candels9/user/khuang/goodss'
vflags_list = {
	0:'normal',
	1:'unresolved',
	2:'could be low-z object, but no proof',
	3:'multiple',
	4:'confirmed a star or low-z interloper',
	10:'diffraction spike or other image artifacts',
	12:'near bright neighbor',
	13:'image quality problems'}

"""
Make cutouts around the specified source ID or position for visual inspection.
Only works for WFC3/IR bands so far (F098M, F105W, F125W, F160W); consider expanding
to the ACS bands?
"""

class cutouts_goodss(object):
	def __init__(self, bands=["f160w"]):
		"""
		Initialize the WCS from headers. bands should be a list (even if there is only one band).
		"""
		self.wcs_dic = {}
		self.bands = bands
		for b in self.bands:
			b = b.lower()
			self.buildwcs(b)
		self.processes = []  # used for storing python.subprocess instances when opening ds9

	def buildwcs(self, band):
		print "Building WCS in %s" % band.upper()
		if band == 'f098m':
			self.image = rootdir+'/gs_all_candels_ers_%s_060mas_v0.5_drz.fits' % band
		else:
			self.image = rootdir+'/gs_all_candels_ers_udf_%s_060mas_v0.5_drz.fits' % band
		hdr = pyfits.getheader(self.image)
		wcs = WCS(hdr)
		self.wcs_dic[band] = wcs

	def id2xy(self, objectid):
		"""
		Converts object ID from the GOODS-S catalog to pixel coordinates. Returned coordinates 
		are in integers.
		objectid is an integer, while band is a string representing the bandpass, e.g., F160W.
		"""
		band = self.bands[0].lower()  # convert ot all lower-case
		ra, dec = c.ra_1[c.id_1==objectid][0], c.dec_1[c.id_1==objectid]
		skycrd = N.array([[ra,dec]])
		pixcrd = self.wcs_dic[band].wcs_sky2pix(skycrd, 1)
		return N.around(pixcrd[0]).astype('int')

	def sky2xy(self, ra, dec, band):
		"""
		Converts (RA, DEC) into pixel coordinates. RA & DEC need to be in degrees.
		"""
		band = band.lower()
		skycrd = N.array([[ra,dec]])
		if band not in wcs_dic.keys():
			buildwcs(band)
		pixcrd = self.wcs_dic[band].wcs_sky2pix(skycrd, 1)
		return N.around(pixcrd[0].astype('int'))

	def xy2temp(self, x, y, band, width=2.0, temp='temp.fits'):
		"""
		For one band only. Make temp.fits
		"""
		#image = rootdir+'/gs_all_candels_ers_udf_%s_060mas_v0.5_drz.fits' % band
		segimage = rootdir+'/gs_all_sx_h_120604_hphotom_comb_seg_psfmatch2h.fits'
		pixwidth = int(N.round(width/0.06))
		xlo = x - pixwidth/2; xhi = x + pixwidth/2
		ylo = y - pixwidth/2; yhi = y + pixwidth/2
		iraf.imcopy(self.image+'[%d:%d,%d:%d]'%(xlo,xhi,ylo,yhi),temp)
		temp_seg = os.path.splitext(temp)
		temp_seg = temp_seg[0]+'_seg.fits'
		iraf.imcopy(segimage+'[%d:%d,%d:%d]'%(xlo,xhi,ylo,yhi),temp_seg)
		return temp, temp_seg

	def xy2temps(self, x, y, width=2.0):
		"""
		For temporary images in self.bands. Calls self.xy2temps for each band.
		"""
		tempfits = ""
		for b in self.bands:
			temp = 'temp_%s.fits' % b
			temp,temp_seg = self.xy2temp(x, y, b, width=width, temp=temp)
			tempfits += '%s ' % temp
			tempfits += '%s ' % temp_seg
		return tempfits

	def cutouts_id(self, objectid, width=2.0):
		"""
		Make square cutouts of width given in arcseconds.
		Specify the object ID from the catalog, and the bands
		one wants the cutouts.
		"""
		if len(glob.glob('temp*.fits')):
			os.system('rm temp*.fits')
		# there are multiple bands; create multiple panels
		x, y = self.id2xy(objectid)
		tempfits = self.xy2temps(x, y, width=width)
		seglow = objectid - 200
		seghigh = objectid + 200
		ds9verb = "ds9 -zscale %s -frame 1 -zoom to fit -frame 2 -scale limits %d %d -zoom to fit " % (tempfits,seglow,seghigh)
		print ds9verb
		#os.system(ds9verb+' &')
		x = subprocess.Popen(ds9verb.split())
		self.processes += [x]

	def clear_processes(self):
		"""
		Clear self.processes
		"""
		for x in self.processes:
			x.kill()


class visual_inspect_goodss(cutouts_goodss):
	"""
	Create cutouts of objects to visually inspect them. The routines that makes the cutouts are
	inherited from the cutouts_goodss class.
	"""
	def __init__(self, id_list, cu, bands=['f160w']):
		self.id = id_list
		self.cu = cu
		self.nobj = len(self.id)
		cutouts_goodss.__init__(self, bands)
		self.vflags = N.ones(len(id_list),'int') * -1
		self.current_index = 0
		self.comments = ['']*len(id_list)
		self.vflags_list = vflags_list
		self.vflags_dic = {'id':self.id, 'vflags':self.vflags}
		self.vflags_dic['current_index'] = self.current_index
		self.vflags_dic['vflags_list'] = self.vflags_list
		self.vflags_dic['comments'] = self.comments
		# Note that vflags_list could be expanded by the user, the the user should register
		# what the additional vflag values mean.

	def load(self, vflagfile):
		"""
		Load the pickled file that contained the dictionary of visual-inspection flags.
		"""
		f = open(vflagfile,'rb')
		self.vflags_dic = cPickle.load(f)
		self.id = self.vflags_dic['id'].copy()
		self.current_index = self.vflags_dic['current_index']
		self.vflags = self.vflags_dic['vflags'].copy()
		self.comments = copy.deepcopy(self.vflags_dic['comments'])
		self.vflags_list = copy.deepcopy(self.vflags_dic['vflags_list'])
		f.close()

	def save(self, vflagfile):
		"""
		Save self.vflags_dic to a pickled file.
		"""
		# first update self.vflags_dic
		self.vflags_dic = {'id':self.id.copy(), 'vflags':self.vflags.copy()}
		self.vflags_dic['current_index'] = self.current_index
		self.vflags_dic['vflags_list'] = vflags_list.copy()
		self.vflags_dic['comments'] = copy.deepcopy(self.comments)
		self.vflags_dic['vflags_list'] = copy.deepcopy(self.vflags_list)
		f = open(vflagfile, 'wb')
		cPickle.dump(self.vflags_dic, f, 2)
		f.close()

	def import_vflags(self,vflagfile):
		"""
		Import visual flags from existing files, to avoid duplicate efforts.
		"""
		f = open(vflagfile, 'rb')
		lines = f.readlines()
		n_updated = 0
		for l in lines:
			x = l.split()
			objid = int(x[0])
			vflag = int(x[1])
			if len(x) > 2:
				comments = ' '.join(x[2:])
			else:
				comments = ''
			#if len(comments) > 0: print comments
			if objid in self.vflags_dic['id']:
				j = N.arange(len(self.vflags_dic['id']))[self.vflags_dic['id']==objid][0]
				#print j
				self.vflags_dic['vflags'][j] = vflag
				self.vflags_dic['comments'][j] = comments
				n_updated += 1

	def n_inspected(self):
		print "%d out of %d objects have been inspected" % (sum(self.vflags_dic['vflags']>=0), len(self.vflags_dic['vflags']))

	def idnow(self):
		"""
		Print the current ID number.
		"""
		indexnow = self.current_index
		return self.id[indexnow]

	def vflagnow(self):
		"""
		Return the vflag of the current ID.
		"""
		indexnow = self.current_index
		return self.vflags[indexnow]

	def vflag_id(self, objectid):
		return self.vflags[self.id==objectid][0]

	def goback(self):
		"""
		Go back to the previous object.
		"""
		self.current_index -= 1
		self.clear_processes()
		print "Current object:"
		print self.idnow()
		print self.vflagnow()

	def next(self):
		"""
		Go to the next object.
		"""
		if self.vflagnow() < 0:
			print "Warning: ID %d still not inspected!!" % self.idnow()
		self.current_index += 1
		self.clear_processes()
		print "Current object:"
		print self.idnow()
		print self.vflagnow()

	def next_uninspected(self):
		"""
		Go to the next un-inspected one (with vflag < 0).
		"""
		while self.vflagnow() >= 0:
			self.next()


	def goto(self, objectid):
		"""
		Go to a specific object with objectid.
		"""
		index_object = N.arange(self.nobj)[self.id==objectid][0]
		self.current_index = index_object
		self.clear_processes()
		print "Current object:"
		print self.idnow()
		print self.vflagnow()

	def recordflag(self, value, comment='', objid=None):
		"""
		Record the visual-inspection flag of objid as value.
		"""
		if objid == None:
			objid = self.idnow()
		indexnow = self.current_index
		self.vflags[indexnow] = value
		self.record_comment(comment)
		if objid != self.idnow():
			pass
		else:
			self.current_index += 1
		if self.current_index >= len(self.id):
			print "Current index is over the length of ID list."
		self.clear_processes()
		print "Current object:"
		print self.idnow()
		print self.vflagnow()

	def record_comment(self, comment_str, objectid=None):
		"""
		Record special comments about the current object.
		"""
		if objectid == None:
			objectid = self.idnow()
		index = N.arange(self.nobj)[self.id==objectid]
		print "index", index
		self.comments[index] = comment_str

	def get_comment(self, objectid):
		"""
		Retrieve the comment of object with ID.
		"""
		index = N.arange(self.nobj)[self.id==objectid]
		print self.comments[index]

	def speczdq_now(self):
		"""
		Print the spec-z quality of this object.
		"""
		dq = self.cu.spec_z_dq[self.cu.id_1==self.idnow()][0]
		print "spec-z quality =", dq

	def SN_now(self,instr='wfc3',band='f160w'):
		"""
		Print the S/N of the given band.
		"""
		sn = (getattr(self.cu,'%s_%s_flux'%(instr,band))/getattr(self.cu,'%s_%s_fluxerr'%(instr,band)))[self.cu.id_1==self.idnow()][0]
		print "S/N in %s %s is %.2f" % (instr,band,sn)

	def inspect_this(self, width=6.0):
		print "Inspecting ID %d" % self.idnow()
		print "spec-z = %.2f" % self.cu.spec_z[self.cu.id_1==self.idnow()][0]
		print "phot-z = %.2f" % self.cu.photo_z[self.cu.id_1==self.idnow()][0]
		self.cutouts_id(self.idnow(), width=width)

	def stat(self):
		n_done = sum(self.vflags >= 0)
		print "%d out of %d objects have been inspected." % (n_done, self.nobj)
		print "A break down of classified objects:"
		flags_dic = {}
		for vf in self.vflags_list.keys():
			if vf not in flags_dic.keys():
				flags_dic[vf] = 0
		for i in range(self.nobj):
			if self.vflags[i] >= 0:
				flags_dic[self.vflags[i]] += 1
		for vf in self.vflags_list.keys():
			print "Flag %d has %d objects." % (vf, flags_dic[vf])

	def writeto(self, catname, comment=False):
		"""
		Write to a text file.
		"""
		f = open(catname,'wb')
		for i in range(self.nobj):
			f.write('%8d %8d  ' % (self.id[i],self.vflags[i]))
			if comment:
				f.write('%s ' % self.comments[i])
			f.write('\n')
		f.close()


