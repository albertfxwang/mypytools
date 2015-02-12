#!/usr/bin/env python

# Plot the specific SFR as a function of redshift
# Use the data compiled from literature:
# Gonzalez et al. 2010
# Smit et al. 2014
# Reddy et al. 2012
import numpy as np
import matplotlib.pyplot as plt
import itertools

scatter_markers = itertools.cycle(['s','v','o','d','p'])
scatter_colors = itertools.cycle(['Red', 'DarkOrange', 'ForestGreen', 'SteelBlue', 'DarkOrchid', 'DarkSlateGray'])


ssfr = {}  # in [z, sSFR, sig(sSFR)] pairs, where sSFR is in Gyr^-1
# sig(sSFR) are in dex; if sig(sSFR) < 0, it is taken as an upper limit
# it's easier to read off the values from Reddy et al. 2012 than from Gonzalez
# et al. 2010...
ssfr['Salim+07'] = [[0.1, 0.15, 0.3]]
ssfr['Noeske+07'] = [[0.3, 0.3, 0.3], [0.55, 0.56, 0.3], [0.75, 0.85, 0.3], [1.0, 1.2, 0.3]]
ssfr['Daddi+07'] = [[2.0, 2.0, 0.32]]
ssfr['Reddy+12'] = [[2.25, 2.4, 0.37], [3.1, 2.3, 0.37]]
ssfr['Bouwens+12'] = [[4.0, 5.1, 0.3], [5.0, 5.0, 0.3], [6.0, 3.1, 0.3], [7.25, 5.3, 0.3]]
# ssfr['Stark+09'] = [[4.0, 2.05, 0.1], [5.0, 2.0, 0.1], [6.0, 2.0, 0.1]]
ssfr['Stark+13'] = [[4.0, 6.0, 0.3], [5.0, 5.9, 0.3], [6.0, 6.2, 0.3], [6.8, 10.0, 0.3]]
ssfr['Labbe+13'] = [[8.0, 11.0, 11.0, 5.0]]
ssfr['Gonzalez+10'] = [[7.2, 2.3, 0.1]]
ssfr['Smit+14'] = [[6.8, 14.0, -0.2]]
# ssfr['Our sample'] is in the format of [zbest, +zerr, -zerr, ssfr_best, +ssfr_err, -ssfr_err]
ssfr_sample = {'MACS1149-JD':[9.4, 0.2, 0.2, 0.88, 18.0, 0.088],
              'MACS1423-1384':[7.4, 0.79, 0.60, 8.7, 12.0, 8.2],
              'RXJ1347-1080':[7.3, 0.28, 0.34, 0.72, 2.3, 0.34],
              'MACS0744-2088':[7.0, 0.3, 0.024, 50., 2.3, 47.0],
              'MACS1423-587':[6.8, 0.05, 5.6, 110., 0.001, 100.],
              'MACS1423-1494':[7.1, 0.2, 0.061, 110, 0.035, 36.],
              'MACS1423-2097':[6.7, 0.15, 0.092, 0.42, 2.4, 0.25],
              'RXJ1347-1216':[6.77, 0.001, 0.001, 110, 0.001, 48.],
              'RXJ1347-1800':[6.7, 0.057, 5.7, 110, 0.001, 100.],
              'RCS2327-1282':[7.1, 0.037, 0.049, 2.8, 0.98, 0.057],
              'MACS0454-1251':[6.1, 0.04, 0.064, 3.4, 22., 0.25],
              'MACS0454-1817':[6.7, 0.0054, 0.098, 2.8, 2.4, 0.7],
              'MACS1149-274':[5.8, 0.09, 0.067, 7.3, 1.2, 6.1],
              'MACS1149-1204':[5.7, 0.17, 0.097, 5.1, 71., 0.92]}
# maybe also plot stellar mass? stellar masses in 10^9 solar mass, corrected for magnification
smass_sample = {'MACS1149-JD': [0.61, 0.19, 0.29],
                'MACS1423-1384': [3.7, 3.5, 2.5],
                'RXJ1347-1080': [1.0, 0.14, 0.65],
                'MACS0744-2088': [1.2, 1.5, 0.091],
                'MACS1423-587':  [0.90, 0.11, 0.49],
                'MACS1423-1494': [0.14, 0.63, 0.0076],
                'MACS1423-2097': [1.9, 1.3, 0.22],
                'RXJ1347-1216':  [0.22, 0.50, 0.03],
                'RXJ1347-1800':  [0.37, 0.82, 0.11],
                'RCS2327-1282':  [2.0, 0.42, 0.46],
                'MACS0454-1251': [4.0, 0.029, 2.2],
                'MACS0454-1817': [12., 2.6, 6.6],
                'MACS1149-274':  [3.3, 1.8, 0.67],
                'MACS1149-1204': [2.8, 0.44, 2.3]}

# sSFR from Dave+11 @ 10^10 M_solar; in [z, log(SSFR)] pairs
ssfr_sim = [[0., -0.9],[1.0, -0.42],[2.2, -0.22],[4., 0.16],[6.0,0.5]]

# marker properties in [marker, facecolor, edgecolor]; default size is 10
markers = {}
up_arrow = [(-3,0),(3,0),(0,0),(0,4),(-1.5,2),(0,4),(1.5,2),(0,4),(0,0)]
markers['Salim+07'] = ['^', 'none', 'RoyalBlue']
markers['Noeske+07'] = ['s', 'none', 'cyan']
markers['Daddi+07'] = ['s', 'blue', 'none']
markers['Reddy+12'] = ['x', 'none', 'green']
markers['Bouwens+12'] = ['o', 'DarkOrange', 'none']
markers['Stark+13'] = ['^', 'magenta', 'none']
markers['Gonzalez+10'] = ['o', 'none', 'purple']
markers['Smit+14'] = [up_arrow, 'none', 'DarkRed']
markers['Labbe+13'] = ['p', 'none', 'CadetBlue']
markers['Our sample'] = ['*', 'black', 'none']
ebar_kwargs = {'mew':2, 'capsize':8, 'ms':10}

def plot_ssfr(filename='/Users/khuang/Dropbox/research/my_papers/surfsup3/figures/ssfr_evolution.eps'):
   # Make a sSFR evolution plot
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_yscale('log')
   i = 0
   # First plot the literature values
   for k in ssfr.keys():
      marker, facecolor, edgecolor = markers[k]
      if facecolor == 'none':
         ecolor = edgecolor
      else:
         ecolor = facecolor
      z_lit = []
      ssfr_val_lit = []
      ssfr_sig_lit = []
      ssfr_low_lit = []
      ssfr_high_lit = []
      for i in range(len(ssfr[k])):
         z_lit += [ssfr[k][i][0]]
         ssfr_val_lit += [ssfr[k][i][1]]
         if len(ssfr[k][i]) == 4:
            # ssfr_sig_lit.append(ssfr[k][i][2:])
            ssfr_high_lit.append(ssfr[k][i][1] + ssfr[k][i][2])
            ssfr_low_lit.append(ssfr[k][i][1] - ssfr[k][i][3])
         else:
            # ssfr_sig_lit += [ssfr[k][i][2]]
            ssfr_high_lit.append(ssfr[k][i][1]*(10.**ssfr[k][i][2]))
            ssfr_low_lit.append(ssfr[k][i][1]*(10.**-ssfr[k][i][2]))
      z_lit = np.array(z_lit)
      ssfr_val_lit = np.array(ssfr_val_lit)
      # ssfr_sig_lit = np.array(ssfr_sig_lit)
      # ssfr_low_lit = ssfr_val_lit * (10.**-ssfr_sig_lit)
      # ssfr_high_lit = ssfr_val_lit * (10.**ssfr_sig_lit)
      if k == 'Smit+14':
         ax.scatter(z_lit, ssfr_val_lit, marker=marker, s=13**2,
                    edgecolor=edgecolor, label=k, linewidths=2)
      else:
         ax.errorbar(z_lit, ssfr_val_lit, 
            yerr=[ssfr_val_lit-ssfr_low_lit,ssfr_high_lit-ssfr_val_lit],
            fmt=marker, mfc=facecolor, mec=edgecolor, mew=2.,
            ecolor=ecolor, label=k, capsize=8, ms=10)

   # Plot simulated sSFR evolution from Dave+11 (vzw model)
   for i in range(len(ssfr_sim)-1):
      if i == 0:
         label = 'Dave+11 (vzw)'
      else:
         label = ''
      ax.plot([ssfr_sim[i][0], ssfr_sim[i+1][0]], 
              [10.**ssfr_sim[i][1], 10.**ssfr_sim[i+1][1]],
              lw=1.0, c='Indigo', label=label)
   # Now plot OUR SAMPLE
   our_marker, our_mfc, our_mec = markers['Our sample']
   set_label = False
   for k in ssfr_sample:
      if not set_label:
         label = 'This work'
         set_label = True
      else:
         label = ""
      z, zerr_hi, zerr_lo = ssfr_sample[k][:3]
      our_ssfr, ssfrerr_hi, ssfrerr_lo = ssfr_sample[k][3:]
      ax.errorbar(z, our_ssfr, 
                  xerr=[[zerr_lo], [zerr_hi]], 
                  yerr=[[ssfrerr_lo], [ssfrerr_hi]],
                  fmt=our_marker, mfc=our_mfc, mec=our_mec,
                  ecolor=our_mfc, label=label, capsize=6, ms=16,
                  elinewidth=1.5, capthick=1.5)

   ax.set_xlabel('Redshift')
   ax.set_ylabel('sSFR ($\mathrm{Gyr}^{-1}$)')
   ax.legend(loc=4, ncol=3, markerscale=0.8, fontsize='large', mode='expand')
   ax.set_xlim(xmin=-0.5, xmax=11)
   ax.set_ylim(ymin=1e-3)
   plt.draw()
   fig.savefig(filename)
   return ax

def plot_ssfr_smass(filename='/Users/khuang/Dropbox/research/my_papers/surfsup3/figures/smass_vs_ssfr.eps'):
   fig = plt.figure()
   ax = fig.add_axes([0.12, 0.15, 0.78, 0.75])
   ax.set_xscale('log')
   ax.set_yscale('log')
   for k in ssfr_sample.keys():
      smass, smass_errhi, smass_errlo = smass_sample[k]
      ssfr, ssfr_errhi, ssfr_errlo = ssfr_sample[k][3:]
      marker = scatter_markers.next()
      color = scatter_colors.next()
      ax.errorbar([smass*1.e9], [ssfr], 
                  xerr=[[smass_errlo*1.e9],[smass_errhi*1.e9]],
                  yerr=[[ssfr_errlo],[ssfr_errhi]], 
                  fmt=marker, mfc=color, mec='none', ecolor=color, 
                  elinewidth=1.2, 
                  capsize=8, capthick=1.2, ms=8, label=k)
   ax.legend(ncol=3, loc='lower center', fontsize='large', mode='expand')
   ax.set_xlabel(r'Stellar mass [$10^9\,M_{\odot}$]')
   ax.set_ylabel(r'sSFR [$\mathrm{Gyr}^{-1}$]')
   ax.set_ylim(ymin=2.e-3, ymax=300)
   ax.set_xlim(xmax=2.e10)
   fig.savefig(filename)
   return ax

