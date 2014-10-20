#!/usr/bin/env python

import transfunc_calc as tfc
from pygoods import *

simcat = '/Users/khuang/Dropbox/Research/bivariate/galfit_transfer/merged_I_091712.fits'
limits = array([[21.0,26.5],[-1.,2.0]])
pixdx = array([0.02,0.02])
cell_widths = array([0.5,0.2])
columns = ['mag','logre']

testkgrid = tfc.kernelgrid(simcat,limits,cell_widths,pixdx,columns,'test.p')
c = testkgrid.simcat
testkgrid.set_detection(c.re_out>0.)
qf = logical_and(((c.re_out_err/c.re_out)<=0.6) , (c.chisqnu<=0.45))
testkgrid.set_inclusion(qf)
testkgrid.set_gtype(c.galaxy_type,disk=0)

testkgrid.transfer_grid_2d()