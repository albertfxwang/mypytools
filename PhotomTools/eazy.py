## Tools for processing/displaying EAZY results
import numpy as np
import os
from pygoods import Ftable, sextractor
from sedtools.bands import filters
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz

eazy_path = '/Users/khuang/bin/eazy'
# a symbol for magnitude upper limits
downarrow = [(-2,0),(2,0),(0,0),(0,-4),(-2,-2),(2,-2),(0,-4),(0,0)]

class PlotEAZY(object):
   def __init__(self, catalog, output_dir='OUTPUT'):
      # assumes that all outputs are within the directory output_dir
      # catalog is the file containing all the fluxes; first column in catalog
      # is the object ID.
      self.c = sextractor(catalog)
      self.objid = self.c._1
      self.output_dir = output_dir
      # also assumes that the first line of catalog is the header line
      f = open(catalog)
      hdr = f.readline()
      f.close()
      self.filter_names = []
      # read filter names; they should be within the keys in 
      # sedtools.bands.filters
      for b in hdr.split():
         if b.lower().startswith('f_'):
            self.filter_names += [b.lower()[2:]]
      # Read in the output photo-z file
      out = sextractor("%s/photz.zout" % self.output_dir)
      hdr = open("%s/photz.zout" % self.output_dir).readline()
      hdr_list = hdr.split()[1:]
      # print hdr_list
      for i in range(1, len(hdr_list)+1):
         setattr(self, hdr_list[i-1].lower(), getattr(out, "_%d" % i))
      self.lambda_max = 30000.  # default value for maximum lambda in plots
      self.lambda_min = 3000.
      self.flux_ujy_max = 1.0

   def plot_Pz(self, objid, ax=None, savefig=True, mode='a', **plot_kwargs):
      assert objid in self.objid, "Object ID %d does not exist in catalog." % objid
      pz = sextractor("%s/%d.pz" % (self.output_dir, objid))
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      prob = np.exp(-0.5 * pz._2)
      probIntegrated = cumtrapz(prob, x=pz._1)[-1]
      prob = prob / probIntegrated
      ax.plot(pz._1, prob, **plot_kwargs)
      ax.set_xlabel('Redshift')
      ax.set_ylabel(r'$P(z)$')
      if mode == 'a':
         z_peak = self.z_a[self.id==objid][0]
      elif mode=='1':
         z_peak = self.z_1[self.id==objid][0]
      ax.set_title('Object %d [z_peak = %.3f]' % (objid, z_peak))
      if savefig:
         fig.savefig("%s/%d_Pz.png" % (self.output_dir, objid))
      return ax, z_peak

   def plot_photom(self, objid, ax=None, savefig=False, ebar_kwargs={'ms':10, 'mec':'blue', 'mfc':'none', 'ecolor':'blue', 'mew':1.1, 'capsize':8, 'capthick':1.5}):
      # Plot the observed fluxes and the best-fit SED
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_yscale('log')
         ax.set_xscale('log')
      obs = sextractor("%s/%d.obs_sed" % (self.output_dir, objid))
      # First, plot the detections
      detect = (obs._2 > 0)
      ax.errorbar(obs._1[detect], obs._2[detect], yerr=obs._3[detect], 
                  fmt='s', **ebar_kwargs)
      ymin = ax.get_ylim()[0]
      # Then plot upper limits
      # yerr_bottom = obs._2[detect==False] - 1.e-10
      # ax.errorbar(obs._1[detect==False], obs._2[detect==False], 
      #             yerr=[yerr_bottom,obs._3[detect==False]], uplims=True, 
      #             **ebar_kwargs)
      ax.scatter(obs._1[detect==False], obs._3[detect==False],
                 marker=downarrow, edgecolor=ebar_kwargs['ecolor'],
                 s=(ebar_kwargs['ms']*2)**2, linewidths=1.2)
      ax.set_ylim(ymin=ymin)
      self.lambda_max = obs._1.max() * 1.2
      self.lambda_min = obs._1.min() * 0.8
      self.flux_ujy_max = obs._2.max() * 10.
      ax.set_xlim(self.lambda_min, self.lambda_max)
      ax.set_ylim(ymax=self.flux_ujy_max)
      return ax

   def plot_SED(self, objid, xmax=None, ax=None, savefig=False, plot_kwargs={'color':'black'}, mode='a'):
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
         ax.set_yscale('log')
         ax.set_xscale('log')
      # Now plot the best-fit template SED
      temp = sextractor("%s/%d.temp_sed" % (self.output_dir, objid))
      l = open("%s/%d.temp_sed" % (self.output_dir, objid)).readlines()[1]
      l = l.split()
      for i in range(len(l)):
         if l[i] == 'z=':
            z_temp = float(l[i+1])
            break
      ax.plot(temp._1, temp._2, **plot_kwargs)
      ax.set_xlim(self.lambda_min, self.lambda_max)
      lambda_ticks = [5000.,10000.,20000.,30000.,40000.]
      ax.set_xticks(lambda_ticks)
      ax.set_xticklabels(map(lambda x: "%.1f" % (x/1.e4), lambda_ticks))
      ax.set_xlabel(r"$\lambda$ [$\mu$m]")
      ax.set_ylabel(r"$F_{\nu}$ [$\mu$Jy]")
      if mode == 'a':
         z_peak = self.z_a[self.id==objid][0]
      elif mode == '1':
         z_peak = self.z_1[self.id==objid][0]
      ax.set_title('Object %d [z_peak = %.3f]' % (objid, z_peak), weight='bold')
      if savefig:
         plt.savefig('%s/%d_SED.png' % (self.output_dir, objid))
      return ax, z_peak

   def plot_all(self, objid, objname="", savefig=True, legend_loc=1, mode='a'):
      fig = plt.figure(figsize=(10,12))
      ax1 = fig.add_subplot(211)
      ax1.set_yscale('log')
      ax1.set_xscale('log')
      ax1 = self.plot_photom(objid, ax=ax1)
      ax1, zp = self.plot_SED(objid, ax=ax1, mode=mode)
      ax2 = fig.add_subplot(212)
      ax2, zp = self.plot_Pz(objid, ax=ax2, savefig=False, mode=mode)
      if len(objname):
         ax1.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
         ax2.set_title('Object %s [z_peak = %.3f]' % (objname, zp))
      plt.draw()
      plt.savefig("%d_SED_Pz.png" % objid)


def plot_HST_IRAC_SED(catalog_hst, catalog_irac, objid, ax=None, colors=['blue', 'black'], savefig=True, legend_loc=2, mode='a'):
   curdir = os.getcwd()
   print "Plotting object %d..." % objid
   os.chdir('hst_only')
   p1 = PlotEAZY(catalog_hst)
   if mode == 'a':
      zp_1 = p1.z_a[p1.id==objid][0]
   elif mode == '1':
      zp_1 = p1.z_1[p1.id==objid][0]
   os.chdir(curdir+'/with_irac')
   p2 = PlotEAZY(catalog_irac)
   if mode == 'a':
      zp_2 = p2.z_a[p2.id==objid][0]
   elif mode == '1':
      zp_2 = p2.z_1[p2.id==objid][0]
   if ax == None:
      fig = plt.figure()
      ax = fig.add_axes([0.12,0.12,0.78,0.78])
   ax.set_yscale('log')
   ax.set_xscale('log')
   ax = p2.plot_photom(objid, ax=ax, savefig=False)
   os.chdir('../hst_only')
   ax, z_peak = p1.plot_SED(objid, ax=ax, savefig=False, mode=mode,
                    plot_kwargs={'color':colors[0],'label':'HST only (z_peak=%.3f)' % zp_1,'lw':1.5}) 
   os.chdir(curdir+'/with_irac')
   ax, z_peak = p2.plot_SED(objid, ax=ax, savefig=False, mode=mode,
                    plot_kwargs={'color':colors[1],'label':'With IRAC (z_peak=%.3f)' % zp_2,'lw':1.5})
   ax.set_xlim(xmin=p2.lambda_min * 0.8, xmax=p2.lambda_max * 1.5)
   ax.legend(loc=legend_loc) 
   os.chdir(curdir)
   if savefig:
      plt.savefig('%d_SED.png' % objid)
   return ax, z_peak

def plot_HST_IRAC_Pz(catalog_hst, catalog_irac, objid, ax=None, colors=['blue', 'black'], savefig=True, legend_loc=1, mode='a'):
   curdir = os.getcwd()
   print "Plotting P(z) for object %d..." % objid
   os.chdir('hst_only')
   p1 = PlotEAZY(catalog_hst)
   if mode == 'a':
      zp_1 = p1.z_a[p1.id==objid][0]
   elif mode == '1':
      zp_1 = p1.z_1[p1.id==objid][0]
   os.chdir(curdir+'/with_irac')
   p2 = PlotEAZY(catalog_irac)
   if mode == 'a':
      zp_2 = p2.z_a[p2.id==objid][0]
   elif mode == '1':
      zp_2 = p2.z_1[p2.id==objid][0]
   if ax == None:
      fig = plt.figure()
      ax = fig.add_axes([0.12,0.12,0.78,0.78])
   os.chdir('../hst_only')
   ax, z_peak = p1.plot_Pz(objid, ax=ax, savefig=False, color=colors[0],
                   label='HST only (z_peak=%.3f)' % zp_1, ls='--', lw=3.0,
                   mode=mode) 
   os.chdir(curdir+'/with_irac')
   ax, z_peak = p2.plot_Pz(objid, ax=ax, savefig=False, color=colors[1],
                   label='With IRAC (z_peak=%.3f)' % zp_2, lw=1.5,
                   mode=mode)
   ax.legend(loc=legend_loc) 
   os.chdir(curdir)
   if savefig:
      plt.savefig('%d_Pz.png' % objid)
   return ax, z_peak

def plot_HST_IRAC_all(catalog_hst, catalog_irac, objid, objname="", colors=['blue','black'], savefig=True, legend_loc=2, SEDtype='default', mode='a'):
   # Plot both P(z) and SED for the object with ID objid
   # catalog_hst: the file name under directory hst_only/
   # catalog_irac: the file name under directory with_irac/
   fig = plt.figure(figsize=(10,12))
   ax1 = fig.add_subplot(2, 1, 1)
   ax1, z_peak = plot_HST_IRAC_SED(catalog_hst, catalog_irac, objid, ax=ax1, 
                           colors=colors, savefig=False, 
                           legend_loc=legend_loc, mode=mode)
   title = 'Object %d (%s) [z_peak = %.3f]' % (objid, objname, z_peak)
   if objname:
      ax1.set_title(title + '\n(%s SED)' % SEDtype)
   ax2 = fig.add_subplot(2, 1, 2)
   ax2, z_peak = plot_HST_IRAC_Pz(catalog_hst, catalog_irac, objid, ax=ax2, 
                          colors=colors, savefig=False, legend_loc=legend_loc,
                          mode=mode)
   if objname:
      ax2.set_title(title)
   if savefig:
      plt.savefig('%d_SED_Pz.png' % objid)
   plt.show()
   return [ax1, ax2]

