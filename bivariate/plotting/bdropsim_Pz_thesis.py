#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import zdist
from my_mplfonts import Menlo, Helvetica

Mbins_default = array([-23.0, -21.0, -19.0])
logRbins_default = array([0.0, 0.4, 0.8])

def binedge(zd):
   """
   Convert lower bin edges in M1500, log10(Re) into strings for easy
   comparison. First round the float numbers to one decimal digit, then
   convert to string.
   """
   M0 = round(zd.M0, 1)
   M0 = str(M0)
   # convert ot string
   logR0 = round(zd.logR0, 1)
   logR0 = str(logR0)
   return M0, logR0

def Pz_bins(zdgrid1, zdgrid2, z0, z1, Mbins=Mbins_default, 
            logRbins=logRbins_default,
            dropname='B-dropouts', figsize=(9,9)):
   """
   Both Mbins and logRbins need to be sorted:
   Mbins: from small (more negative) to large.
   logRbins: from small to large.
   Plot zdgrid1 for both GOODS (blue solid) and zdgrid2 for UDF (red dashed) 
   depths.
   """
   # I used zdgrid/zdgrid_bdrops_goods_130524.p and 
   # zdgrid/zdgrid_bdrops_udf_130524.p for zdgrid2
   Mbins = array(Mbins)
   Mbins = around(Mbins, 1)
   logRbins = array(logRbins)
   logRbins = around(logRbins, 1) 
   x1 = searchsorted(zdgrid1.Mbins, Mbins) - 1
   y1 = searchsorted(zdgrid1.logrbins, logRbins) - 1
   x2 = searchsorted(zdgrid2.Mbins, Mbins) - 1
   y2 = searchsorted(zdgrid2.logrbins, logRbins) - 1
   fig = plt.figure(figsize=figsize)
   nx = len(Mbins)
   ny = len(logRbins)
   print "nx, ny", nx, ny
   grid = AxesGrid(fig, (0.1,0.1,0.8,0.8), 
                   nrows_ncols=(ny,nx), 
                   axes_pad=0.02,
                   share_all=False,
                   aspect=False,
                   label_mode="L")
   #axes = {}
   for i in range(len(Mbins)):
      for j in range(len(logRbins)):
         # reverse the order in logRbins
         k = i + ny * j
         #ax = fig.add_subplot(nx, ny, k)
         #axes[k] = ax
         ax = grid[k]
         print k
         M, logR = ("%.1f" % Mbins[i], "%.1f" % logRbins[len(logRbins)-j-1])
         #for zk in zdgrid1.zdist.keys():
         #Mk, logRk = binedge(zdgrid1.zdist[zk])
         #if (Mk, logRk) == (M, logR):
         #if zk == (x1[i], y1[j]):
         zk1 = (x1[i], y1[j])
         zk2 = (x2[i], y2[j])
         zd1 = zdgrid1.zdist[zk1]
         zd2 = zdgrid2.zdist[zk2]
         # the zdist instance for this bin
         ax.plot(zd1.zarr, zd1.Pz, '-', c='blue', lw=2.0)
         ax.plot(zd2.zarr, zd2.Pz, '--', c='red', lw=2.0)
         ax.set_xlim(z0, z1)
         ax.set_ylim(0., 1.0)
         ax.set_xticks(arange(z0, z1, 0.5))
         ax.set_yticks(arange(0.2, 1.0, 0.2))
         if k < nx:
            # To show range of M1500
            ax.set_title('$%.1f\leq M_{1500}<%.1f$' % \
                        (float(M),float(M)+0.5),
                        size=14)
         if (k+1) % nx == 0:
            # the right-most axes: show logR bin range
            ax.text(1.08, 0.5, '$%.1f\leq \logR_e < %.1f$' % \
                    (float(logR),float(logR)+0.2),
                    transform=ax.transAxes,
                    ha='center', va='center',
                    size=14, rotation='vertical')
            #ax.yaxis.set_label_position('right')
            #ax.set_ylabel('$%.1f\leq \logR_e < %.1f$' % \
            #              (float(logR),float(logR)+0.2),
            #              size=14)
            #ax.set_ylabel('Y label', size=14)
               
         if k % nx == 0:
            ax.set_ylabel('P(z)', font_properties=Helvetica(18))
            yticks = ax.get_yticks()
            ax.set_yticklabels(map(lambda x: "%.1f" % x, yticks),
                               font_properties=Helvetica(14))
         #else:
         #   ax.set_yticklabels('')

         if k > (nx * (ny-1) - 1):
            xticks = ax.get_xticks()
            ax.set_xticklabels(xticks, font_properties=Helvetica(14))
            ax.set_xlabel('Redshift', font_properties=Helvetica(18))
         #if k == (nx * ny -1):
         #   ax.set_xticklabels(arange(3.0, 5.5, 0.5), 
         #                      font_properties=Menlo(12))
         #else:
         #   ax.set_xticklabels('')

         continue

   fig.suptitle('$P(z)$ for %s' % dropname, font_properties=Helvetica(22))
