#!/usr/bin/env python

import numpy as np
# from pyraf import iraf
import os, sys
# iraf.stsdas()
# iraf.stsdas.analysis.dither()
import pywcs, pyfits
import scipy
from scipy.ndimage import zoom

def cutout_30mas_v1(h_seg, v_drz):
   """
   Because the v1.0, 30mas version of the F606W mosaic includes the parallel 
   fields and therefore covers a larger area than the CANDELS footprint 
   (defined by the v0.5 mosaic), I am making a cutout from the 30mas mosaic
   to cover the same area as the v0.5 60mas mosaics.
   """
   hdr1 = pyfits.getheader(h_seg)
   hdr2 = pyfits.getheader(v_drz)
   nx1 = hdr1['naxis1']
   ny1 = hdr1['naxis2']
   # Now calculate the corners of the cutout in the 30mas frame; 1=60mas frame,
   # 2 = 30mas frame
   wcs1 = pywcs.WCS(hdr1)
   wcs2 = pywcs.WCS(hdr2)
   sky00 = wcs1.wcs_pix2sky([[1, 1]], 1)
   corner00 = np.floor(wcs2.wcs_sky2pix(sky00, 1)).astype('int')[0]
   sky11 = wcs1.wcs_pix2sky([[nx1, ny1]], 1)
   corner11 = np.ceil(wcs2.wcs_sky2pix(sky11, 1)).astype('int')[0]
   xlo, ylo = corner00
   xhi, yhi = corner11
   print "xlo, xhi, ylo, yhi", xlo, xhi, ylo, yhi
   output = os.path.splitext(v_drz)[0] + '_center.fits'
   v_drz_array = pyfits.getdata(v_drz)
   v_drz_hdr = pyfits.getheader(v_drz)
   v_drz_hdr['crpix1'] = v_drz_hdr['crpix1'] - xlo
   v_drz_hdr['crpix2'] = v_drz_hdr['crpix2'] - ylo
   v_drz_array_new = v_drz_array[ylo:yhi+1, xlo:xhi+1]
   pyfits.append(output, v_drz_array_new, v_drz_hdr)
   # iraf.imcopy("%s[%d:%d,%d:%d]" % (v_drz, xlo, xhi, ylo, yhi), output)


def segmap_to_30mas(segimage):
   """
   Drizzle the 60mas seg map onto the 30mas pixel grid.
   whtimage is the weight of the segmentation map, and should be 1.0
   everywhere.
   """
   output = os.path.splitext(segimage)[0] + '_30mas.fits'
   if os.path.exists(output):
      os.remove(output)
   hdr = pyfits.getheader(segimage)
   seg = pyfits.getdata(segimage)
   nx = hdr['NAXIS1']
   ny = hdr['NAXIS2']
   # iraf.wdrizzle(data=segimage, outdata=output, outweig='weight_seg_60mas.fits',
   #               in_mask=whtimage, outnx=2*nx, outny=2*ny, kernel='point',
   #               scale=0.5, xsh=0.0, ysh=0.0, shft_un='output', 
   #               shft_fr='output')
   # iraf.imlintran(input=segimage, output=output, xrotation=0., yrotation=0.,
   #                xmag=0.5, ymag=0.5, ncols=2*nx, nlines=2*ny, 
   #                interpolant='drizzle[1.0]')
   seg2 = zoom(seg, 2.0, order=0).astype('int16')
   print seg2.shape
   # Now update the WCS header keywords
   # wcs = pywcs.WCS(hdr)
   hdr['crpix1'] = hdr['crpix1'] * 2
   hdr['crpix2'] = hdr['crpix2'] * 2
   hdr['cd1_1'] = hdr['cd1_1'] / 2.
   hdr['cd2_2'] = hdr['cd2_2'] / 2.
   pyfits.append(output, seg2, hdr)
   return output

def match_segmap_goods_v2_footprint(segimage, drzimage_v2, output=None):
   """
   Match the footprint of the CANDELS segmentation map to that of the GOODS/ACS
   v2 mosaics.
   """
   hdr1 = pyfits.getheader(segimage)
   hdr2 = pyfits.getheader(drzimage_v2)
   wcs1 = pywcs.WCS(hdr1)
   wcs2 = pywcs.WCS(hdr2)
   sky00 = wcs2.wcs_pix2sky([[1, 1]], 1)
   corner00 = wcs1.wcs_sky2pix(sky00, 1)[0]
   corner00 = np.around(corner00).astype('int')
   nx2 = hdr2['naxis1']
   ny2 = hdr2['naxis2']
   sky11 = wcs2.wcs_pix2sky([[nx2, ny2]], 1)
   corner11 = wcs1.wcs_sky2pix(sky11, 1)[0]
   corner11 = np.around(corner11).astype('int')
   xlo1, ylo1 = corner00
   xhi1, yhi1 = corner11
   print "[xlo:xhi, ylo:yhi]", xlo1, xhi1, ylo1, yhi1
   newshape = (xhi1-xlo1+1, yhi1-ylo1+1)
   assert newshape == (nx2, ny2), "Shape of new seg map does not match the shape of the input drz image..."
   # Now make a cutout of the seg map
   if output==None:
      output = os.path.splitext(segimage)[0] + '_m2goods.fits'
   if os.path.exists(output):
      os.remove(output)
   seg = pyfits.getdata(segimage)
   seg2 = seg[ylo1:yhi1+1, xlo1:xhi1+1]
   # hdr1['crpix1'] = hdr1['crpix1'] - xlo1
   # hdr1['crpix2'] = hdr1['crpix2'] - ylo1
   hdr1['crpix1'] = hdr2['crpix1']
   hdr1['crpix2'] = hdr2['crpix2']
   print "shape of GOODS/ACS v2 mosaic: [%d, %d]" % (hdr2['naxis1'], hdr2['naxis2'])
   print "new shape of the segmap: [%d, %d]" % (seg2.shape[1], seg2.shape[0])
   pyfits.append(output, seg2, hdr1)


if __name__ == "__main__":
   segimage = sys.argv[1]
   print "segimage = ", segimage
   segmap_to_30mas(segimage)
   print "Done."


