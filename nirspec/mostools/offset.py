"""
Determine the offset between multi-slit mask images.
"""

import spectools,special_functions,scipy
import pyfits

from scipy import stats

def findoffset(scidata,coord,offset,slice=None):
	"""
	findoffset(scidata,coord,offset,slice=None)

	Finds the center of a starbox.

	Inputs:
	  scidata - input starbox slit
	  coord   - An array or polynomial describing straightening transform
	  offset  - y-offset for the straightening
	  slice   - sub-slice to use

	Outputs:
	  center of star in starbox
	"""

	"""
	First we straighten the starbox, then find the brightest pixel and
	  choose an aperture (in the x-direction) around that point. The median
	  of this aperture is our star profile that we center wrt.
	"""
	straight = spectools.resampley(scidata,coord,offset,slice=slice)
	illum = stats.stats.median(straight,axis=0)
	peak = illum.argmax()
	start = peak-100
	end = peak+100
	if start<0:
		start = 0
	if end>illum.size:
		end = illum.size
	slice = stats.stats.median(straight[:,start:end],axis=1)
	xvals = scipy.arange(slice.size)*1.

	""" Sharpen the stellar profile. """
	slice = slice**4

	sum = slice.sum()
	slice *= xvals
	center = slice.sum()/sum
	return center
