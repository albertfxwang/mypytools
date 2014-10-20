#!/usr/bin/env python

from numpy import *
from dropsim_util import plot_completeness as pc
import matplotlib.pyplot as plt

arcsec_z2 = 8.597640108504383  # physical distance in kpc at z=2 corresponding to 1"

c1 = pc.plt_comp('udropsim_goodss_udf_130417.fits')
c2 = pc.plt_comp('udropsim_goodss_deep_130417.fits')
c3 = pc.plt_comp('udropsim_goodss_wide_130417.fits')

#fig = plt.figure()
#ax = fig.add_subplot(111)

re_in_lim_kpc = array([0.3,0.5,1.0,1.5,2.0,2.5,3.0])  # in kpc
re_in_lim_as = re_in_lim_kpc/arcsec_z2  # in arcsec
re_in_lim = re_in_lim_as/0.06  # in pixels
logre_in_lim = log10(re_in_lim)
mag_in_array = arange(21.,31.5,0.5)
# Only plot devaucouleurs profile---gtype=1
ocrit1 = c1.gtype==1
ocrit2 = c2.gtype==1
ocrit3 = c3.gtype==1
# Plot 1 --- all size bins in UDF
lines = []
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
l1 = c1.plot_comp_fixed_re('H',ax=ax1,logre_in_lim=[-2.,logre_in_lim[0]],
   label='Re<0.3 kpc',othercrit=ocrit1,pltkw={'lw':1.2})
lines += [l1]
for i in range(1,len(logre_in_lim)-1):
   logre0 = logre_in_lim[i]
   logre1 = logre_in_lim[i+1]
   label = "%.1f kpc < Re < %.1f kpc" % (re_in_lim_kpc[i],re_in_lim_kpc[i+1])
   li = c1.plot_comp_fixed_re('H',ax=ax1,logre_in_lim=[logre0,logre1],label=label,
      mag_in_array=mag_in_array,othercrit=ocrit1,pltkw={'lw':1.2})
ax1.legend(loc=3)
plt.title('UDF depth',size=16)
# Plot 2 --- all size bins in GOODS-S Deep
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
l1 = c2.plot_comp_fixed_re('H',ax=ax2,logre_in_lim=[-2.,logre_in_lim[0]],
   label='Re<0.3 kpc',othercrit=ocrit2,pltkw={'lw':1.2})
lines += [l1]
for i in range(1,len(logre_in_lim)-1):
   logre0 = logre_in_lim[i]
   logre1 = logre_in_lim[i+1]
   label = "%.1f kpc < Re < %.1f kpc" % (re_in_lim_kpc[i],re_in_lim_kpc[i+1])
   li = c2.plot_comp_fixed_re('H',ax=ax2,logre_in_lim=[logre0,logre1],label=label,
      mag_in_array=mag_in_array,othercrit=ocrit2,pltkw={'lw':1.2})
ax2.legend(loc=3)
plt.title('GOODS-S Deep',size=16)
# Plot 3 --- all size bins in GOODS-S Wide
fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
l1 = c3.plot_comp_fixed_re('H',ax=ax2,logre_in_lim=[-2.,logre_in_lim[0]],
   label='Re<0.3 kpc',othercrit=ocrit3,pltkw={'lw':1.2})
lines += [l1]
for i in range(1,len(logre_in_lim)-1):
   logre0 = logre_in_lim[i]
   logre1 = logre_in_lim[i+1]
   label = "%.1f kpc < Re < %.1f kpc" % (re_in_lim_kpc[i],re_in_lim_kpc[i+1])
   li = c3.plot_comp_fixed_re('H',ax=ax3,logre_in_lim=[logre0,logre1],label=label,
      mag_in_array=mag_in_array,othercrit=ocrit3,pltkw={'lw':1.2})
ax3.legend(loc=3)
plt.title('GOODS-S Wide',size=16)
# Plot 4 --- selected size bins in all depths
fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
logrelim = logre_in_lim[1:3]; print logrelim
line11 = c1.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label='0.5 kpc < Re < 1.0 kpc',
   mag_in_array=mag_in_array,pltkw={'ls':'-','color':'blue','lw':1.2},othercrit=ocrit1)
line12 = c2.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label="",
   mag_in_array=mag_in_array,pltkw={'ls':'-.','color':'blue','lw':1.2},othercrit=ocrit2)
line13 = c3.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label="",
   mag_in_array=mag_in_array,pltkw={'ls':':','color':'blue','lw':1.2},othercrit=ocrit3)  
logrelim = logre_in_lim[3:5]; print logrelim
line21 = c1.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label='1.5 kpc < Re < 2.0 kpc',
   mag_in_array=mag_in_array,pltkw={'ls':'-','color':'green','lw':1.2},othercrit=ocrit1)
line22 = c2.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label="",
   mag_in_array=mag_in_array,pltkw={'ls':'-.','color':'green','lw':1.2},othercrit=ocrit2)
line23 = c3.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label="",
   mag_in_array=mag_in_array,pltkw={'ls':':','color':'green','lw':1.2},othercrit=ocrit3)
logrelim = logre_in_lim[5:]; print logrelim
line31 = c1.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label='2.5 kpc < Re < 3.0 kpc',
   mag_in_array=mag_in_array,pltkw={'ls':'-','color':'red','lw':1.2},othercrit=ocrit1)
line32 = c2.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label="",
   mag_in_array=mag_in_array,pltkw={'ls':'-.','color':'red','lw':1.2},othercrit=ocrit2)
line33 = c3.plot_comp_fixed_re('H',ax=ax4,logre_in_lim=logrelim,label="",
   mag_in_array=mag_in_array,pltkw={'ls':':','color':'red','lw':1.2},othercrit=ocrit3)
ax4.legend(loc=3)
ax4.set_title('GOODS-S SExtractor simulations',size=16)
ax4.text(0.05,0.85,'solid:UDF\ndot-dashed:Deep\ndotted:Wide',transform=ax4.transAxes,fontsize=12,
   verticalalignment='top',bbox=dict(facecolor='wheat',alpha=0.5))
