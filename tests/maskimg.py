#!/usr/bin/env python

from numpy import *
import pyfits
import os, glob

# try creating masks using robust statistics
# more specifically, do not use sigma-clipping; instead use estimates
# of scale such as MAD or MADN (normalized median absolute deviation)

def mask_madn(sci_img, flg_img, out_img, f=3.0, maxflag=64):
	# sci_img: input image name to be masked
	# flg_img: the flag image
	# maxflag: threshold of flag value, above which we mask the image
	h1 = pyfits.open(sci_img)
	imgshape = shape(h1[0].data)
	img1 = h1[0].data.ravel()
	h2 = pyfits.open(flg_img)
	img2 = h2[0].data.ravel()
	if not len(glob.glob(out_img)):
		os.system('cp %s %s' % (sci_img, out_img))
		#os.remove(out_img)
	h3 = pyfits.open(out_img, mode='update')
	# now do the masking
	flag = (img2<=maxflag)
	median1 = median(compress(flag,img1))
	MADN = median(abs(compress(flag,img1) - median1))/0.675  # calculates MADN
	print "MADN = ", MADN
	#mask = zeros(len(img1))
	# mask where pixel deviates from median by f MADN or where flag > maxflag
	mask = where((abs(img1-median1)>=(f*MADN))|(img2>maxflag), 1, 0)
	mask = reshape(mask, imgshape)
	#img3 = pyfits.PrimaryHDU(mask)
	#h3 = pyfits.HDUList([img3])
	#h3.writeto(out_img)
	h3[0].data = mask
	h3.flush()
	h3.close()
	h1.close()
	h2.close()
	print "Done."
	
	
	