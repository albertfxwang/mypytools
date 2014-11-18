#!/usr/bin/env python

# Plot the specific SFR as a function of redshift
# Use the data compiled from literature:
# Gonzalez et al. 2010
# Smit et al. 2014
# Reddy et al. 2012
import numpy as np
import matplotlib.pyplot as plt

ssfr = {}  # in [z, sSFR, sig(sSFR)] pairs, where sSFR is in Gyr^-1
# sig(sSFR) are in dex; if sig(sSFR) < 0, it is taken as an upper limit
# it's easier to read off the values from Reddy et al. 2012 than from Gonzalez
# et al. 2010...
ssfr['Salim+07'] = [[0.1, 0.15, 0.3]]
ssfr['Noeske+07'] = [[0.3, 0.3, 0.3], [0.55, 0.56, 0.3], [0.75, 0.85, 0.3], [1.0, 1.2, 0.3]]
ssfr['Daddi+07'] = [[2.0, 2.0, 0.32]]
ssfr['Reddy+12'] = [[2.25, 2.4, 0.37], [3.1, 2.3, 0.37]]
ssfr['Bouwens+11'] = [[4.0, 5.1, 0.3], [5.0, 5.0, 0.3], [6.0, 3.1, 0.3], [7.25, 5.3, 0.3]]
ssfr['Stark+09'] = [[4.0, 2.05, 0.1], [5.0, 2.0, 0.1], [6.0, 2.0, 0.1]]
ssfr['Gonzalez+10'] = [[7.2, 2.3, 0.1]]
ssfr['Smit+14'] = [[6.8, 14.0, -0.2]]
ssfr['Our sample'] = [[9.37, 3.17, 0.3], [6.02, 4.84, 0.3], [6.01, 4.35, 0.3],
                      [6.10, 2.57, 0.3], [7.64, 2.31, 0.3], [5.90, 5.38, 0.3],
                      [8.57, 2.57, 0.3], [8.4, 2.57, 0.3], [7.61, 38.08, 0.3],
                      [7.22, 3.53, 0.3], [6.61, 3.17, 0.3], [6.20, 1.87, 0.3],
                      [6.48, 1.87, 0.3], [6.75, 2.31, 0.3], [6.03, 8.19, 0.3],
                      [7.14, 2.08, 0.3], [6.74, 2.08, 0.3], [7.16, 2.08, 0.3]]
# marker properties in [marker, facecolor, edgecolor]; default size is 10
markers = {}
up_arrow = [(-3,0),(3,0),(0,0),(0,4),(-1.5,2),(0,4),(1.5,2),(0,4),(0,0)]
markers['Salim+07'] = ['^', 'none', 'RoyalBlue']
markers['Noeske+07'] = ['s', 'none', 'cyan']
markers['Daddi+07'] = ['s', 'blue', 'none']
markers['Reddy+12'] = ['x', 'none', 'green']
markers['Bouwens+11'] = ['o', 'DarkOrange', 'none']
markers['Stark+09'] = ['^', 'magenta', 'none']
markers['Gonzalez+10'] = ['o', 'none', 'purple']
markers['Smit+14'] = [up_arrow, 'none', 'DarkRed']
markers['Our sample'] = ['o', '0.5', 'none']
ebar_kwargs = {'mew':2, 'capsize':8, 'ms':10}

def plot_ssfr(our_sample=True):
   # Make a sSFR evolution plot
   fig = plt.figure()
   ax = fig.add_subplot(111)
   ax.set_yscale('log')
   i = 0
   for k in ssfr.keys():
      if (k == 'Our sample') and (our_sample==False):
         continue
      marker, facecolor, edgecolor = markers[k]
      if facecolor == 'none':
         ecolor = edgecolor
      else:
         ecolor = facecolor
      z = []
      ssfr_val = []
      ssfr_sig = []
      for i in range(len(ssfr[k])):
         z += [ssfr[k][i][0]]
         ssfr_val += [ssfr[k][i][1]]
         ssfr_sig += [ssfr[k][i][2]]
      z = np.array(z)
      ssfr_val = np.array(ssfr_val)
      ssfr_sig = np.array(ssfr_sig)
      for i in range(len(ssfr[k])):
         if i == 0:
            label = k
         else: 
            label = ""
         if ssfr_sig[i] > 0:
            ssfr_low = ssfr_val[i] * (10.**-ssfr_sig[i])
            ssfr_high = ssfr_val[i] * (10.**ssfr_sig[i])
            if k == 'Our sample':
               ax.errorbar(z[i], [ssfr_val[i]], 
                        yerr=[[ssfr_val[i]-ssfr_low],[ssfr_high-ssfr_val[i]]],
                        fmt=marker, mfc=facecolor, mec=edgecolor, 
                        ecolor=ecolor, label=label, capsize=5, ms=6)
            else:
               ax.errorbar(z[i], [ssfr_val[i]], 
                        yerr=[[ssfr_val[i]-ssfr_low],[ssfr_high-ssfr_val[i]]],
                        fmt=marker, mfc=facecolor, mec=edgecolor, 
                        ecolor=ecolor, label=label, **ebar_kwargs)
         else:
            ax.scatter(z, ssfr_val, marker=marker, facecolor=facecolor, 
                       edgecolor=edgecolor, label=label, s=18**2, 
                       linewidths=2.0)


   ax.set_xlabel('z')
   ax.set_ylabel('sSFR (Gyr$^{-1}$)')
   ax.legend(loc=4, ncol=2, markerscale=0.8)
   ax.set_xlim(xmin=-0.5)
   plt.draw()
   return ax
