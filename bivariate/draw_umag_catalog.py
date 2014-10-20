#!/usr/bin/env python

from numpy import *
import udropsim_uflux as uu
from pygoods import Ftable
import fitsutil
import os

rootdir = '/Users/khuang/Dropbox/Research/bivariate/udrops_fitting'
rootcat = rootdir+'udropsim_run2m_goodss_deep_130505.fits'
c = Ftable(rootcat)

def draw_umags(niter=1):
	for i in range(niter):
		vimos_u_mag = uu.draw_umag(c.u_mag_in,
                 'simcatalogs/udrops/uvimos_mag_pdf.p')
		vimos_u_mag_1sig = uu.draw_umag(c.u_mag_in,'simcatalogs/udrops/uvimos_mag_1sigma_pdf.p',onesigma=True)
		vimos_u_mag = minimum(vimos_u_mag,vimos_u_mag_1sig)
		os.system('cp %s udropsim_umag_drawn_%d.fits'%(rootcat,i))
		fitsutil.add_columns('udropsim_umag_drawn_%d.fits'%i, ['vimos_u_mag_out'],[vimos_u_mag],['D'])
		
