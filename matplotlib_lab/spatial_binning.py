#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from numpy import random

# Example taken from the following page:
# http://www.physics.ucdavis.edu/~dwittman/Matplotlib-examples/

# Example 1: hexbin, from http://matplotlib.org/examples/pylab_examples/hexbin_demo.html
np.random.seed(0)
n = 100000
x = np.random.standard_normal(n)
y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
xmin = x.min()
xmax = x.max()
ymin = y.min()
ymax = y.max()

plt.subplots_adjust(hspace=0.5)
plt.subplot(121)
plt.hexbin(x,y, cmap=plt.cm.YlOrRd_r)
plt.axis([xmin, xmax, ymin, ymax])
plt.title("Hexagon binning")
cb = plt.colorbar()
cb.set_label('counts')

plt.subplot(122)
plt.hexbin(x,y,bins='log', cmap=plt.cm.YlOrRd_r)
plt.axis([xmin, xmax, ymin, ymax])
plt.title("With a log color scale")
cb = plt.colorbar()
cb.set_label('log10(N)')

plt.show()

# Example 2: square bins, using numpy.histogram2d
hist,xedges,yedges = numpy.histogram2d(x,y,bins=40,range=[[-6,4],[-4,6]])
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
imshow(hist.T,extent=extent,interpolation='nearest',origin='lower')
colorbar()
show()