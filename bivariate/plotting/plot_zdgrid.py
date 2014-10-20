#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import zdist
import cPickle, time


def flash_sel_comp(zdgrid, z0, z1, dz, field='UDF', t=1.0):
   """
   Flash the selection completeness in redshift intervals, with time steps
   spefified by t.
   """
   fig = plt.figure()
   ax = fig.add_subplot(111)
   cbar = 0
   for z in arange(z0, z1, dz):
      zdgrid.calc_sel_comp(z0=z, z1=z+dz, show=False)
      cax = ax.imshow(zdgrid.sel_comp_map.swapaxes(0,1),
         extent=(0,shape(zdgrid.sel_comp_map)[0],0,shape(zdgrid.sel_comp_map)[1]),
         vmin=0., vmax=1.0)
      if cbar == 0:
         plt.colorbar(cax)
         cbar = 1
      ax.set_title('z between [%.1f, %.1f] in %s' % (z, z+dz, field), size=18)
      plt.draw()
      time.sleep(t)

