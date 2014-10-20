from nirspec import niredux as nir
# To reduce our NIRSPEC on on 2014/09/01 on MACS2129
# Last updated: 14/09/09
# For the low resolution mode of NIRSPEC, pixel scale is 0.144" in the spatial
# direction and 0.193" in the spectral direction.

# KH.371 is ~ 9.57" away from the trace along the slit --> ~ 66.5 pixels away
datadir = '/Users/khuang/Dropbox/Research/Keck/NIRSPEC/2014sep02_B/spec'
flats = [datadir+'/sep02s0006.fits',
         datadir+'/sep02s0007.fits',
         datadir+'/sep02s0008.fits',
         datadir+'/sep02s0009.fits',
         datadir+'/sep02s0010.fits']

# Use the trace in A370 slit as the star exposure, after aggressively masking
# out the bright sky lines. It seems to work somewhat?
# Here is how I made the mask:
# In [1]: exp29 = pyfits.open('sep02s0029.fits',ignore_missing_end=True)[0].data
# In [2]: exp29mask = np.where(exp29>100, 1, 0)
# 
# In [3]: exp29masked = np.where(exp29mask>0, np.median(exp29), exp29)
# 
# In [4]: pyfits.append('sep02s0029masked.fits',exp29masked)
# stars = [datadir+'/sep02s0029masked.fits']
stars = ['/Users/khuang/Dropbox/Research/Keck/NIRSPEC/Raw/aug14s0019.fits']  
# use Chris F's trace image

imgpairs = [[datadir+'/sep02s0001.fits', datadir+'/sep02s0002.fits'],
            [datadir+'/sep02s0003.fits', datadir+'/sep02s0004.fits']]

nir.niredux(flats,stars,imgpairs,'macs2129',userlow=9470.,userhigh=11210.)

