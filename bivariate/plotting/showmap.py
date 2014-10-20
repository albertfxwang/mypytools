#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt

cmap_hot = mpl.cm.hot

def showmap(marray, limits, bw=array([0.5,0.2]), 
   dtype='float', vrange=None, mask=None, aspect=1.0,
   newfig=True, fig=None, ax=None, cmap=cmap_hot, cbar_extend='neither'):
   # marray -- the 2-D array to be shown
   # limits -- the lower and upper limits of the (x, y) axes; should be [[xlo,xhi],[ylo,yhi]]
   # bw -- bin width in (dx, dy)
   # dtype -- data type of marray
   # if provide a mask, mask should have the same shape as marray
   # vrange -- the dynamic range used to show marray (default is automatically determined)
   # aspect -- aspect ratio of the axes (y/x)
   # newfig -- to create a new figure window or not (default = yes)
   # fig -- the figure in which we show marray IF newfig == False
   # ax -- the axis instance in which we show marray IF newfig == False
   # cmap -- matplotlib color map
   # cbar_extend -- how to format the color bar for out-of-dynamic-range values
   nbins = (limits[:,1] - limits[:,0]) / bw
   nbins = nbins.astype('int')
   if not (shape(marray) == nbins).all():
      print "shape of input array is [%d, %d], incompatible with number of bins [%d, %d]" % (
         shape(marray)[0], shape(marray)[1], nbins[0], nbins[1])
      raise ValueError
   xbins = arange(limits[0][0], limits[0][1], bw[0])
   ybins = arange(limits[1][0], limits[1][1], bw[1])
   xticklabels = []; yticklabels = []
   for k in range(len(xbins)):
      xticklabels += ['%.1f' % xbins[k]]
   for k in range(len(ybins)):
      yticklabels += ['%.1f' % ybins[k]]
   if mask != None:
      masked = zeros(npix, 'bool')
            
   if newfig == True:
      fig = plt.figure()
      ax = fig.add_subplot(111)
   if vrange==None:
      vrange=[min(marray.ravel()), max(marray.ravel())]
   if mask != None:
      marray = ma.masked_where(masked==True, marray)
   cax = ax.imshow(marray.swapaxes(0,1), origin='lower', vmin=vrange[0], aspect=aspect,
      vmax=vrange[1], cmap=cmap, interpolation='none')
   cbar = fig.colorbar(cax, ax=ax, extend=cbar_extend)
   ax.set_xticks(arange(0,nbins[0])); ax.set_xticklabels(xticklabels)
   ax.set_yticks(arange(0,nbins[1])); ax.set_yticklabels(yticklabels)
   
   return marray, fig, ax, cbar