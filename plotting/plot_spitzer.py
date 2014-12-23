## From Marusa (2014/12/18)

#from setup_matplotlib import *
import numpy as np
#from setup_matplotlib import *
import matplotlib.pyplot as plt
from matplotlib import rcParams, rc
rc('text.latex', preamble='\usepackage{sfmath}')
rc('text.latex', preamble='\usepackage{amsmath}')
rc('text', usetex=True)
#####Note, ony SWARPed images work (WCS!!)
import aplpy
import alpha
import Image
import sys
import pywcs
import pyfits
import numpy as np

filters = ['RC_drz']
ra = [333.73852,333.73925,333.73682]
dec = [-14.003327,-13.999389,-14.003195]
delta = [400/3600.0, 30/3600.0, 10/3600.0]
for index in range(len(filters)):
        fits = filters[index]+'.fits'
        jpeg = filters[index]+'.pdf'
        for c in range(len(ra)):
            
# Show the RGB image
            if (index == 0): 
                fig = aplpy.FITSFigure(fits)
                fig.show_rgb('MACS2214_ds9.jpg')
            else:
                fig = aplpy.FITSFigure(fits)
                fig.show_grayscale(vmin = 0, vmax = scales[index], stretch='arcsinh')
 
# Overlay a grid
#fig.add_grid()
#fig.grid.set_alpha(0.1)

# Save image
            jpeg = filters[index]+'_'+str(c)+'.pdf'

            fig.recenter(ra[c], dec[c], width=delta[c], height=delta[c])
            #fig.show_markers(ra + delta/np.cos(dec)/50, dec, layer='marker_set_1', edgecolor='red',
            #                         facecolor='none', marker=1, s=2000, alpha=1,lw=4)
            #fig.show_markers(ra, dec+delta/20, layer='marker_set_2', edgecolor='red',
            #                         facecolor='none', marker=2, s=2000, alpha=1, lw=4)
            fig.axis_labels.set_font(size='xx-large', weight='medium', \
                         stretch='normal', family='sans-serif', \
                         style='normal', variant='normal')
                         #fig.add_label(ra - 0.2*delta/np.cos(dec)/2,dec+1.08*delta/2,filters[index], size=36, fontweight='bold', color='black')
            #fig.add_grid()
            #fig.show_regions('irfov.reg')
            #fig.grid.set_alpha(0.5)
            fig.ticks.set_color('black')
            # Save image
            fig.set_theme('publication')
            if (c > 0):
                fig.axis_labels.hide()
                fig.tick_labels.hide()
            fig.save(jpeg, dpi=100)
            fig.close()
        
