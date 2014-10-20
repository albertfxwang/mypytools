"""
Determine the spectral resolution.
"""

import scipy
from scipy import ndimage

import special_functions

def clipped(data,clip=3.5):
        d = data.copy()
        mean,std = d.mean(),d.std()

        while 1:
                size = d.size
                d = d[abs(d-mean)<clip*std]
                if d.size==size or d.size==0:
                        return mean,std
                mean,std = d.mean(),d.std()
        return mean,std


def measure(spectrum,num=25):
	"""
	measure(spectrum,num=25)

	Deterines the spectral resolution for a spectrum by measuring the width
	  of arc or skylines.

	Inputs:
	  spectrum - 1d spectrum
	  num      - maximum number of lines to measure (not implemented)

	Outputs:
	  median resolution of all lines extracted
	"""

	data = spectrum.astype(scipy.float64)
	size = data.size

	pixels = scipy.arange(size)*1.

	avg,std = clipped(data,3.)
	thresh = avg + 30.*std  # Only measure 30 sigma lines

	""" Identify peaks. """
	mask = ndimage.maximum_filter(data,13)
	mask = scipy.where(mask==data,data,0)
	lines = scipy.where((mask>thresh)&(mask<50000.))[0]
	count = 0

	vals = []
	for pos in lines:
		""" Don't use lines near the edges. """
		if pos<5 or pos+6>size:
			continue
		fitdata = scipy.zeros((11,2))
		fitdata[:,0] = pixels[pos-5:pos+6]
		fitdata[:,1] = data[pos-5:pos+6]

		par = scipy.zeros(4)
		par[1] = data[pos]
		par[2] = pixels[pos]
		par[3] = 1.

		fit,chi2 = special_functions.ngaussfit(fitdata,par)
		fit[2] -= fitdata[0,0]
		"""
		Reject fits that were 'too' wide or narrow, or not near the
		  expected center.
		"""
		if fit[3]>4. or fit[3]<0.8 or fit[2]<2 or fit[2]>9:
			continue

		vals.append(fit[3])

		count += 1
	vals = scipy.asarray(vals)
	if vals.size==0:
		return 1.
	return scipy.median(vals)
