#!/usr/bin/env python

from numpy import *
from pygoods import *
import KPDFadaptnumpy as KPDF
import draw_ndist as dn
import matplotlib as mpl
import matplotlib.pyplot as plt
import cPickle
import mlutil
from scipy.interpolate import RectBivariateSpline

"""
A general-purpose script that calculates transfer-function kernels in cells 
of the parameter space, using data from simulation output.
Many routines are inherited from adaptivetransfer2d.py. For now only supports 
2D kernels.
"""

### Change Logs ###
# 13/06/06: fixed the bug in transfer_grid_2d that didn't use outindex to
#				select objects for kernel construction (lines 286-287).

__author__ = "Kuang-Han Huang"


class kernel(object):
	""" A convolution kernel associated with a certain range of parameter space
	"""
	def __init__(self,x0_cell,x1_cell,pixdx,completeness,kernel_array):
		self.x0_cell = x0_cell                         
		# Where to start using this kernel
		self.x1_cell = x1_cell                         
		# Where to stop using this kernel
		self.xc_cell = (x0_cell + x1_cell) / 2.
		# The coordinate at the center of this cell
		self.pixdx = array(pixdx)                         
		# kernel pixel size (a length-2 array)
		self.i0 = -array(shape(kernel_array))/2    
		# the pixel index of the lower-left corner, relative to center
		self.i1 = array(shape(kernel_array))/2     
		# the pixel index of the upper-right corner, relative to center
		self.rx0 = -array(shape(kernel_array))*pixdx/2.  
		# the pixel coordinate of the lower-left corner, relative to center
		self.rx1 = array(shape(kernel_array))*pixdx/2.   
		# the pixel coordinate of the upper-right corner, relative to center
		self.completeness = completeness   
		# completeness is specified separately from kernel_array
		self.kernel = kernel_array   
		# always normalized to sum(kernel_array.ravel()) = 1
		
	def update(self,kernel_array,pixdx=None):
		if pixdx == None:
			pass
		else:
			self.pixdx = pixdx
		self.i0 = -array(shape(kernel_array))/2
		self.i1 = array(shape(kernel_array))/2
		self.rx0 = -array(shape(kernel_array))*self.pixdx/2.
		self.rx1 = array(shape(kernel_array))*self.pixdx/2.
		#self.completeness = sum(kernel_array.ravel()) * self.pixdx.prod()
		self.kernel = kernel_array
		
	def show(self,imshow_kw={'origin':'lower'}):
		plt.imshow(self.kernel.swapaxes(0,1),**imshow_kw)
		plt.title(self.index)
		
	def __getstate__(self):
		odict = self.__dict__.copy()
		return odict
	
	def __setstate__(self,d):
		self.__dict__.update(d)

   
class kernelgrid(object):
	""" 
	An evenly spaced grid of convolution kernels 
	For each kernel grid, one needs the following attributes:
	-- source simulation catalog
	-- two columns in the catalog that form the 2-D parameter space
	   -- will look for columns with names XXX_in and XXX_out in the catalog, where XXX
	      are the column names specified here
	-- an array that specifies the inclusion of points (used for calculating completeness)
	-- limits of the entire parameter space; a 2x2 array
	-- cell widths (used to calculate the number of cells), a length-2 array
	-- "pixel" size of the distribution
	-- kernel shape of the kernel PDF estimator (Gaussian or Epanechnikov)
	"""
	def __init__(self,simcat,limits,cell_widths,pixdx,filename,
	             columns=['mag','re'],kpdf_shape="gaussian",
	             xmax=array([2.56,2.56])):
		"""
		simcat:			file name of the GALFIT simulation catalog
		limits: 			a 2x2 array of the limits in input (mag, log10(Re)) for the 
				  			kernel grid
		cell_widths: 	a 1x2 array for the bin widths of each kernel
		pixdx: 			a 1x2 array for the **pixel** sizes of the kernel
		columns: 		the attribute names of the two columns (default to 
		               ['mag','re'])
		filename: 		the output kernel file name
		kpdf_shape: 	the kernel shape (Gaussian or Epanechnikov) used in kernel 
		            	PDF estimation
		xmax: 			the extent of each kernel
		"""
		self.simcatname = simcat
		self.limits=limits
		ncells = (limits[:,1]-limits[:,0])/cell_widths
		if (around(ncells) % 1 == [0.0,0.0]).all():
			self.ncells = ncells.astype('int')
		else:
			self.ncells = ncells.astype('int') + 1 # call the applied region of each kernel "cell"
		self.cell_widths = cell_widths  # cell size
		self.pixdx = pixdx   # pixel size
		self.kernels = {}
		self.kpdf_shape = kpdf_shape  # shape of kernel used in kernel PDF estimation
		self.xmax = xmax  # the extent (in parameter-space coordinate) of each kernel
		self.filename = filename
		self.columns = columns
		# an array to store discrete sampling of completeness at each cell center
		self.completeness_array = zeros(self.ncells)
		# the x (magnitude) coordinate at the center of each cell
		self.xc_array = arange(limits[0][0]+cell_widths[0]/2.,
		                       limits[0][1], cell_widths[0])
		# the y (logRe) coordinate at the center of each cell
		self.yc_array = arange(limits[1][0]+cell_widths[1]/2.,
		                       limits[1][1], cell_widths[1])
		# read the simulation catalog
		if self.simcatname.endswith('.fits'):
			self.simcat = Ftable(self.simcatname)
			self.nobj_simcat = len(self.simcat.d)
			self.simcat_fmt = 'fits'
			self.col_in_1 = self.simcat.__getitem__('%s_in'%self.columns[0])
			self.col_out_1 = self.simcat.__getitem__('%s_out'%self.columns[0])
			self.col_in_2 = self.simcat.__getitem__('%s_in'%self.columns[1])
			self.col_out_2 = self.simcat.__getitem__('%s_out'%self.columns[1])
			self.include = ones(len(self.simcat.d),'bool')  # for completeness
			# self.include will be set with another method; if not set then assume all points
			# in the catalog are included (100% completeness)
		else:
			self.simcat = sextractor(simcat)
			self.nobj_simcat = len(self.simcat)
			self.simcat_fmt = 'sex'
			self.col_in_1 = self.simcat.__getattribute__('%s_in'%self.columns[0])
			self.col_out_1 = self.simcat.__getattribute__('%s_out'%self.columns[0])
			self.col_in_2 = self.simcat.__getattribute__('%s_in'%self.columns[1])
			self.col_out_2 = self.simcat.__getattribute__('%s_out'%self.columns[1])
			self.include = ones(len(self.simcat),'bool')  # for completeness
		self.gtype = zeros(self.nobj_simcat, 'int')  # a flag for galaxy type...default is that disk = 0
		self.disk = 0   # the value in self.gtype that represents disks
		self.include = ones(self.nobj_simcat,'bool')
		self.detect = ones(self.nobj_simcat,'bool')

		
	
	def __getnewargs__(self):
		return (self.simcatname,self.limits,self.cell_widths,self.pixdx,
		        self.columns,self.filename,self.kpdf_shape)
			
	def set_detection(self, det_array):
		if len(det_array) != self.nobj_simcat:
			raise ValueError, "the length of det_array is different from the \
				number of sources in the catalog."
		
		self.detect = det_array
		self.nobj_detect = sum(det_array)
		
	def set_inclusion(self, inc_array):
		if len(inc_array) != self.nobj_simcat:
			raise ValueError, "the length of inc_array is different from the \
				number of sources in the catalog."
		
		a = (inc_array==True)
		b = (self.detect==True)
		self.include = logical_and(a,b)  # fold in detection
		self.nobj_include = sum(a)
		
	def set_gtype(self, gtype_array, disk=0):
		self.gtype = gtype_array
		self.disk = disk
		
	def set_other_attributes(self,attr_dic):
		for k in attr_dic.keys():
			setattr(self,k,attr_dic[k])
	
	def getkernel(self,x):  
		""" 
		Return the kernel corresponding to a given set of parameter-space 
		coordinates.
		x -- a LIST of parameter-space coordinates
		"""
		y = (x-self.limits[:,0]) / self.cell_widths
		i = floor(y).astype("int")
		if self.kernels.has_key(tuple(i)):
			return self.kernels[tuple(i)]
		else:
			print "No kernel for ",x
			return -99

	def transfer_grid_2d(self,ntot=-1,diskfrac=0.5,other_attr={},verbose=0,
	                     expand=[0.,0.]):
		#simcat,limits,binwd,pixdx,xmax,verbose=1,
		#kernel_estimator="gaussian",chinu_max=0.4,reerrrmax=0.6,
		#diskfrac=0.7,ntot=60000,sersicnbins=arange(0.,9.),method='galfit'):
		"""
		This method takes simulation input & output and construct 2-D transfer 
		function kernel grid.
		Supposed to be a GENERAL PURPOSE script, and the parameter-space 
		coordinates are specified by self.col_[in/out]_[12].
		ntot: 		total number of simulation points drawn to make kernels
		diskfrac: 	the fraction of disks in the objects that are used to 
						calculate the kernels
		other_attr: other attributes one would like to tack onto the kernelgrid 
						instance (mostly for documentation purposes)
		
		Notes:
		completeness in each cell is calculated by 
		sum(self.include) / len(self.include) IN EACH CELL. One needs to call 
		the set_inclusion method to set up self.include.
		""" 
		#nbins = around((limits[:,1]-limits[:,0])/binwd).astype('int')
		print "Kernel estimator form:", self.kpdf_shape
		print "Kernel grid limits:", self.limits
		print "Number of kernel:", self.ncells
		print "Kernel width:", self.xmax
		#if ntot==None: ntot = self.nobj_include
		

		last_invcovariance = ones(shape(self.limits),"float")

		# Draw points from simulation to match Sersic n distribution of LBGs
		# Then clean up using GALFIT quality flags
		c = self.simcat
		self.expand = expand
		#gtype = c.galaxy_type
		#nout = c.sersicnout
		if ntot < 0:
			outindex = arange(len(self.simcat.d))
		else:
			outindex = dn.draw_ndist(diskfrac,ntot,self.gtype,disk=self.disk)
		print "len(outindex)", len(outindex)
		siminput = array([self.col_in_1.take(outindex),
		                  self.col_in_2.take(outindex)])
		print "shape(siminput)",shape(siminput)
		simoutput = array([self.col_out_1.take(outindex),
		                   self.col_out_2.take(outindex)])
		#residuals = siminput - simoutput
		residuals = simoutput - siminput

		## INITIALIZE KERNEL GRID
		if len(other_attr)>0:
			self.set_other_attributes(other_attr)

		kmaxbins = around(self.xmax/self.pixdx).astype('int')

		# Create an array of bin coordinates in the kernel grid
		fx0,fx1,fxi = mlutil.bincoords2d(self.limits,self.cell_widths)   
		# x0,x1,xi will be Nx2 arrays, where N is the total number of objects
		if verbose > 1:
			print "fx0 = ",fx0
			print "fx1 = ",fx1
			print "fxi = ",fxi
 
		# Create list of bin ranges and indices in the array of kernels
		# fx0,fx1,fxi will be 2-d arrays
		# fx0[i] will be the coordinates of the lower-left point of the ith grid
		# fx1[i] will be the coord of the upper-right point of the ith grid
		# fxi[i] will be the index of the ith grid
		
		# Iterate through each bin, calculate completeness, 
		# clean up using GALFIT quality flags, then make kernels
		for i in range(shape(fx0)[0]):
			# setup bin boundaries and indices
			bx0_cell = fx0[i].copy()   
			# Lower bound for data used to make this kernel
			bx1_cell = fx1[i].copy()   
			# Upper bound
			bx0_sw = fx0[i].copy()  
			# Copies if bin boundaries are expanded
			bx1_sw = fx1[i].copy()
			bxi = fxi[i]
			# Expand the cells a little to use more data... to create kernels
			# that vary more smoothly across the grid.
			bx0_cell = bx0_cell - array(expand)
			bx1_cell = bx1_cell + array(expand)
			print "bx0_cell", bx0_cell
			print "bx1_cell", bx1_cell

			#binx0_cell = outer(bx0_cell,ones(self.nobj_simcat))
			#binx1_cell = outer(bx1_cell,ones(self.nobj_simcat))
			binx0_cell = outer(bx0_cell, ones(len(outindex)))
			binx1_cell = outer(bx1_cell, ones(len(outindex)))

			# determines which input points are in the given bin
			withincell = (siminput >= binx0_cell) & (siminput < binx1_cell)  
			# get data within bin
			withincell = multiply.reduce(withincell)  
			# multiply every element in withinbin together
			print "shape(withincell)",shape(withincell)
			n_withincell = sum(withincell)
			print "sum(withincell)", n_withincell

			# CALCULATE COMPLETENESS
			# completeness = (num. of objects detected by SE & pass GALFIT 
			#                  quality check) / (total num. of objects detected 
			#						 in that cell)
			cell_detect = withincell & self.detect.take(outindex)
			cell_include = withincell & self.include.take(outindex)
			n_detected = sum(cell_detect)
			n_included = sum(cell_include)  
			# number of sources included in this bin
			if n_included > 4:
				completeness = float(n_included) / float(n_detected)
			else:
				completeness = 0.0
			print "completeness,",completeness

			input_cell_1 = compress(cell_include,siminput[0,:])                   
			input_cell_2 = compress(cell_include,siminput[1,:])
			res_cell_1 = compress(cell_include,residuals[0,:])
			res_cell_2 = compress(cell_include,residuals[1,:])
			input_cell = array([input_cell_1,input_cell_2])
			if len(input_cell_1) > 0:
				print input_cell_1.min(), input_cell_1.max()
				print input_cell_2.min(), input_cell_2.max()
			res_cell = array([res_cell_1,res_cell_2])
			nmodel_cell = n_included

			# INITIALIZE KERNEL PROPERTIES
			karray = zeros(kmaxbins,"float")           # Create an empty kernel
			center = self.xmax/2. # the index of center of kernel
			# e.g. if kernel is 128 by 128 array, then index of center will be [64,64]
			kx = arange(-self.xmax[0]/2.+self.pixdx[0]/2.,self.xmax[0]/2.+self.pixdx[0]/2.,self.pixdx[0])
			ky = arange(-self.xmax[1]/2.+self.pixdx[1]/2.,self.xmax[1]/2.+self.pixdx[1]/2.,self.pixdx[1])
			gdata = KPDF.MPDF2DGrid2Array(kx,ky) # output PDF grid
			print "shape(gdata)",shape(gdata)
			covariance=cov(res_cell,rowvar=1)

			if n_included <= 4:  # less than 4 points in this bin
				print "*** WARNING -- kernel has less than 4 points***"
				print "bx0_bin, bx1_bin = ",bx0_cell,bx1_cell
				pdf2d = ones((shape(kx)[0],shape(ky)[0]),"float")
				print "shape(pdf2d)", shape(pdf2d)
				lam = None
				covariance = None
				bandwidth = -1.
			else: 
				rr = transpose(res_cell)
				if verbose > 1:
					print "kx=",kx
				#covariance = cov(rr)

				# Compute scale factors to use for KPDF 
				# This is a kludge to prevent crashing on singular covariance matrices
				try:  
					invcovariance = linalg.inv(covariance)
				except: 
					print "*** Warning, singular inverse convariance matrix, trying Moore-Penrose inverse ***" 
					try:
						invcovariance = linalg.pinv(covariance,rcond=1.e-10)
					except:
						# This will bomb if it happens on the first pass
						print "*** No joy...using last good invcovariance ***" 
						invcovariance = last_invcovariance
				last_invcovariance = invcovariance

				renorm_constant = sqrt(linalg.det(covariance))
				bandwidth=KPDF.MPDFOptimumBandwidth(rr)

				if (self.kpdf_shape=="epanechnikov"):
					lam = ones((shape(gdata)[0]),"float")
					pdf2d = KPDF.MPDFEpanechnikov(rr,rr,bandwidth,lam,invcovariance,renorm_constant)
					lam = KPDF.MPDFAdapt(pdf2d)
					pdf2d = KPDF.MPDFEpanechnikov(rr,gdata,bandwidth,lam,
						invcovariance,renorm_constant)
					pdf2d = reshape(pdf2d,(shape(kx)[0],shape(ky)[0]))
				else:
					lam = ones((shape(gdata)[0]),"float")
					pdf2d = KPDF.MPDFGaussian(rr,rr,bandwidth,lam,
						invcovariance,renorm_constant)
					lam = KPDF.MPDFAdapt(pdf2d)
					pdf2d = KPDF.MPDFGaussian(rr,gdata,bandwidth,lam,
						invcovariance,renorm_constant)
					pdf2d = reshape(pdf2d,(shape(kx)[0],shape(ky)[0]))

			if verbose > 1:
				print "pdf2d=",pdf2d

			if sum(pdf2d.ravel()) > 0.:
				print "sum(pdf2d.ravel())",sum(pdf2d.ravel())
				#karray = normalize(pdf2d,1.0)
				karray = pdf2d / sum(pdf2d)
				#karray = pdf2d
				print "sum(karray.ravel())", sum(karray.ravel())
			else:
				karray = zeros(shape(pdf2d))

			# Normalize the kernel (accounting for incompleteness)
			if n_included <= 4: # if number of points used in kernel is less than 4
				completeness = 0.
			karray = karray * completeness
			print "sum(karray.ravel())", sum(karray.ravel())
			if verbose:
				for i in range(len(bx0_cell)):
					print "%7.3g %7.3g |" % (bx0_cell[i],bx1_cell[i]),
				print "%7d %7.3f " % (n_included,completeness)
				print '\n'

			# Recover cell boundaries for book-keeping
			bx0_cell = bx0_cell + expand
			bx1_cell = bx1_cell - expand
			kclass = kernel(bx0_cell,bx1_cell,self.pixdx,completeness,karray)
			kclass.ninput_cell = n_withincell
			kclass.nmodel_cell = n_included # number of points recovered by SE in the bin
			kclass.ndetect_cell = n_detected
			kclass.covariance = covariance
			kclass.lam = lam
			kclass.bandwidth = bandwidth
			kclass.index = tuple(bxi)
			kclass.withincell = withincell
			#kclass.res_cell = res_cell
			# Add the kernel to the grid
			self.kernels[tuple(bxi)] = kclass
			# Set self.completeness_array
			self.completeness_array[bxi[0],bxi[1]] = completeness
		#delattr(self, 'simcat')
		print "len(self.kernels)", len(self.kernels) 
	
	def interpolate_completeness(self, limits, pixdx):
		"""
		Given two arrays denoting the x and y coordinates of each pixel, 
		interpolate the GALFIT completeness at each pixel.
		"""
		# limits should be in the format of (x0, x1, y0, y1)
		# make sure that limits is within self.limits
		xpix = arange(limits[0][0], limits[0][1], pixdx[0])
		ypix = arange(limits[1][0], limits[1][1], pixdx[1])
		limits = limits.ravel()
		self.comp_interp = RectBivariateSpline(self.xc_array, self.yc_array,
		                                       self.completeness_array,
		                                       bbox=self.limits.ravel())
		cm = self.comp_interp(xpix, ypix)
		cm = minimum(cm, 1.0)
		cm = maximum(cm, 0.0)
		self.completeness_model = cm


	def __getstate__(self):
		# defines how this class will be pickled. Necessary to pickle protocol 2?
		odict = self.__dict__.copy()
		if 'simcat' in odict.keys():
			del odict['simcat']
		return odict
	
	#def __setstate__(self, d):
		#if d['simcatname'].endswith('.fits'):
		#	d['simcat'] = Ftable(d['simcatname'])
		#else:
		#	d['simcat'] = sextractor(d['simcatname'])
		#self.__dict__.update(d)
		
	#def Save(self):
	#	if self.filename==None:
	#		raise ValueError, "Please define the file name of this kernel grid."
	#	f = open(self.filename,'wb')
	##	f.close()
		
	#def Load(self,filename=None):
	#	if filename==None:
	#		filename = self.filename
	#	f = open(filename, 'rb')
	#	tmp_dict = cPickle.load(f)
	#	f.close()
	#	self.__dict__.update(tmp_dict)
		
		
