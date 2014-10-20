#!/usr/bin/env python

from numpy import *
import bivRL as bl
import matplotlib as mpl
import matplotlib.pyplot as plt
from my_mplfonts import Helvetica

p1 = array([-1.7, -21.0, 0.7, 0.7, 0.3])
limits = array([[25.5, 27.0], [0.2, 1.0]])
pixdx = array([0.02, 0.02])
mci = bl.mconvert('/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/kcorr/M1500_to_f775w.txt')

def bivmodel_pix_bin(param=p1, limits=limits, pixdx=pixdx):
   """
   Plot a schematic diagram demonstrating the concept of "pixels" and "bins"
   for the bivariate size-luminosity model distribution.
   """
   model = bl.bivariate_RL(param, limits, pixdx, 'f435w', 'goods',
                           kgrid=None, zdgrid=None, mc=mci,
                           add_interloper=False)

   xmax, ymax = shape(model.model)
   fig1 = plt.figure()
   cax = plt.imshow(model.model.swapaxes(0,1), origin='lower', 
              extent=(0,xmax,0,ymax), interpolation='none')
   cax.axes.set_xticks(arange(0, xmax))
   cax.axes.set_yticks(arange(0, ymax))
   xpix = around(arange(limits[0][0], limits[0][1], pixdx[0]),2)
   ypix = around(arange(limits[1][0], limits[1][1], pixdx[1]),2)[1:]
   xticks = map(lambda x:'%.2f'%x, arange(limits[0][0],limits[0][1],0.5))
   yticks = map(lambda y:'%.2f'%y, arange(limits[1][0],limits[1][1],0.2))
   #print xticks
   xticklabels = map(lambda x:'%.1f'%x if '%.2f'%x in xticks else '',
                     xpix)
   xticklabels[-1] = '%.1f' % (limits[0][1])
   yticklabels = map(lambda y:'%.1f'%y if '%.2f'%y in yticks else '',
                     ypix)
   yticklabels[-1] = '%.1f' % (limits[1][1])
   cax.axes.set_xticklabels(xticklabels, font_properties=Helvetica(14))
   cax.axes.set_yticklabels(yticklabels, font_properties=Helvetica(14))
   plt.xlabel('magnitude', font_properties=Helvetica(18))
   plt.ylabel('log10(Re) [pixel]', font_properties=Helvetica(18))
   # Now draw a white box around a bin
   x0bin = int(round(0.5/pixdx[0])); x1bin = int(round(1.0/pixdx[0]))
   y0bin = int(round(0.2/pixdx[1])); y1bin = int(round(0.4/pixdx[1]))
   #cax.axes.plot([x0bin, x1bin], [y0bin, y0bin], '-', lw=4, color='white')
   cax.axes.fill([x0bin, x1bin, x1bin, x0bin],
                 [y0bin, y0bin, y1bin, y1bin],
                 ec='white', fill=True, fc='white', lw=4, ls='solid',
                 alpha=0.8)
   cax.axes.text((x0bin+x1bin)/2., y0bin-0.1/pixdx[1], 
                 '$\Delta m=0.5$ mag', font_properties=Helvetica(20),
                 horizontalalignment='center', verticalalignment='bottom',
                 color='black')
   cax.axes.text((x0bin-0.1/pixdx[0]), (y0bin+y1bin)/2., 
                 '$\Delta\log R_e=0.2$', font_properties=Helvetica(20),
                 rotation='vertical', color='black',
                 horizontalalignment='left', verticalalignment='center')
   cax.axes.text((x0bin+x1bin)/2., (y0bin+y1bin)/2., 
                 'bin', font_properties=Helvetica(22),
                 color='black', horizontalalignment='center',
                 verticalalignment='center')
   # Draw a white pixel
   x0pix = 1.0 / pixdx[0]; x1pix = x0pix + 1
   y0pix = 0.6 / pixdx[1]; y1pix = y0pix + 1
   cax.axes.fill([x0pix,x1pix,x1pix,x0pix], [y0pix,y0pix,y1pix,y1pix],
                 fc='white', ec='white', fill=True)
   cax.axes.annotate('pixel', xy=[x0pix+0.5,y0pix+0.5], 
                     xytext=[x0pix+5,y0pix-5],
                     arrowprops=dict(arrowstyle="->",
                     connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                     font_properties=Helvetica(16))

   plt.grid(b=True, which='major', color='0.5', linewidth=0.5, linestyle='-')
   plt.grid(b=True, which='minor', color='0.5', linewidth=0.5, linestyle='-')
