#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

def scatter_hist(x, y, binwidth, xlabel='', ylabel='', scatter_label='', legend_loc=0, figsize=(8,8), **scatter_kw):
   """
   Make a scatter plot with histograms on each axis.
   Taken from the example at 
   http://matplotlib.org/examples/pylab_examples/scatter_hist.html
   """
   from matplotlib.ticker import NullFormatter
   nullfmt   = NullFormatter()         # no labels
   fig = plt.figure(figsize=figsize)
   # enforce x and y to be numpy arrays
   x = np.array(x)
   y = np.array(y)
   
   # definitions for the axes
   left, width = 0.1, 0.65
   bottom, height = 0.1, 0.65
   bottom_h = left_h = left+width+0.02
   
   rect_scatter = [left, bottom, width, height]
   rect_histx = [left, bottom_h, width, 0.2]
   rect_histy = [left_h, bottom, 0.2, height]
      
   axScatter = plt.axes(rect_scatter)
   axHistx = plt.axes(rect_histx)
   axHisty = plt.axes(rect_histy)
   
   # no labels
   axHistx.xaxis.set_major_formatter(nullfmt)
   axHisty.yaxis.set_major_formatter(nullfmt)
   
   # the scatter plot:
   axScatter.scatter(x, y, label=scatter_label, **scatter_kw)
   axScatter.legend(loc=legend_loc, scatterpoints=1)
   
   # now determine nice limits by hand:
   if type(binwidth) == type(1.0):
      xbinwidth = binwidth
      ybinwidth = binwidth
   else:
      xbinwidth, ybinwidth = binwidth
   xbuffer = (x.max() - x.min()) / 10.
   ybuffer = (y.max() - y.min()) / 10.
   axScatter.set_xlim( x.min() - xbuffer, x.max() + xbuffer )
   axScatter.set_ylim( y.min() - ybuffer, y.max() + ybuffer )
   
   xbins = np.arange(x.min(), x.max() + xbinwidth, xbinwidth)
   hx = axHistx.hist(x, bins=xbins, histtype='step', lw=2)
   yticks = axHistx.get_yticks()
   yticks = yticks[1:-1][::2]
   axHistx.set_yticks(yticks)
   ybins = np.arange(y.min(), y.max() + ybinwidth, ybinwidth)
   hy = axHisty.hist(y, bins=ybins, orientation='horizontal', histtype='step',
                     lw=2)
   xticks = axHisty.get_xticks()
   xticks = xticks[1:-1][::2]
   axHisty.set_xticks(xticks)
   # axHisty.set_xticklabels(['']+map(lambda x:'%d'%int(round(x)), xticks[1:-1])+[''], 
   #                         rotation=45)
   
   axHistx.set_xlim( axScatter.get_xlim() )
   axHisty.set_ylim( axScatter.get_ylim() )

   if len(xlabel):
      axScatter.set_xlabel(xlabel)
   if len(ylabel):
      axScatter.set_ylabel(ylabel)
   
   plt.show()

   return axScatter, axHistx, axHisty