#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import plot_comp
from pygoods import Ftable
from my_mplfonts import Helvetica, Menlo

"""
Calculate the detection completeness map in (F850LP, logRe) plane.
"""
#df = myfonts.Menlo

rootdir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit'
#c1 = Ftable(rootdir+'/simcatalogs/bvdrops/bvdropsim_run8_nolya_111011.fits')
# GOODS depth
#c2 = Ftable(rootdir+'/simcatalogs/bvdrops/bvdropsim_udf_nolya_out.fits')
# HUDF depth
#acsz_v2_dir = '/Users/khuang/Dropbox/Research/catalogs/goods/north/acsz_v2'
#udf_v2_dir = '/Users/khuang/Dropbox/Research/catalogs/udf'
#cz_goods = Ftable(acsz_v2_dir+'/h_goods_nz_r2.0z.fits')
#cz_udf = Ftable(udf_v2_dir+'/h_udf_z.fits')

def make_compmap(c1, c2):
   # Now calculate the detection completenesses
   compmap1 = plot_comp.sexcompmap(c1, [20.,30.], [-0.6,1.8], 0.5, 0.2)
   compmap2 = plot_comp.sexcompmap(c2, [20.,30.], [-0.6,1.8], 0.5, 0.2)

   dim = shape(compmap1)
   # Plot
   fig = plt.figure(figsize=(8,8))
   ax1 = fig.add_subplot(211)
   # Show GOODS-Depth map first
   cax1 = ax1.imshow(compmap1.swapaxes(0,1), origin='lower', vmin=0., vmax=1.,
                     extent=[0,dim[0],0,dim[1]])
   ax1.set_xticks(arange(0,dim[0],2)[1:])
   ax1.set_yticks(arange(0,dim[1],2))
   ax1.set_xticklabels(arange(20.,30.5,1.)[1:], font_properties=Helvetica(14))
   ax1.set_yticklabels(arange(-0.6,2.0,0.4), font_properties=Helvetica(14))
   ax1.set_ylabel('log10(Re) [pixels]', font_properties=Helvetica(18))
   ax1.text(0.5, 0.1, 'GOODS-Depth', font_properties=Helvetica(18), color='white',
            transform=ax1.transAxes, verticalalignment='center',
            horizontalalignment='center')
   cb1 = fig.colorbar(cax1)
   cb1.set_ticks(arange(0., 1.2, 0.2))
   ax1.set_xlabel('input magnitude', font_properties=Helvetica(18))

   # Show HUDF-Depth map 
   ax2 = fig.add_subplot(212)
   cax2 = ax2.imshow(compmap2.swapaxes(0,1), origin='lower', vmin=0., vmax=1.,
                     extent=[0,dim[0],0,dim[1]])
   ax2.set_xticks(arange(0,dim[0],2)[1:])
   ax2.set_yticks(arange(0,dim[1],2))
   ax2.set_xticklabels(arange(20.,30.5,1.)[1:], font_properties=Helvetica(14))
   ax2.set_yticklabels(arange(-0.6,2.0,0.4), font_properties=Helvetica(14))
   ax2.set_ylabel('log10(Re) [pixels]', font_properties=Helvetica(18))
   ax2.text(0.5, 0.1, 'UDF-Depth', font_properties=Helvetica(18), color='white',
            transform=ax2.transAxes, verticalalignment='center',
            horizontalalignment='center')
   ax2.set_xlabel('input magnitude', font_properties=Helvetica(18))
   cb2 = fig.colorbar(cax2)
   cb2.set_ticks(arange(0., 1.2, 0.2))

   fig.suptitle('Detection Completeness in F850LP Mosaics',
                font_properties=Helvetica(24))
