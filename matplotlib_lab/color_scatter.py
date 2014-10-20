#!/usr/bin/env python
# Example taken from http://stackoverflow.com/questions/10761429/how-to-modify-2d-scatterplot-to-display-color-based-off-third-array-in-csv-file

import matplotlib.pyplot as plt
from matplotlib  import cm
import numpy as np

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
ax.set_title("X vs Y AVG",fontsize=14)
ax.set_xlabel("XAVG",fontsize=12)
ax.set_ylabel("YAVG",fontsize=12)
ax.grid(True,linestyle='-',color='0.75')
x = np.random.random(30)
y = np.random.random(30)
z = np.random.random(30)

# scatter with colormap mapping to z value
ax.scatter(x,y,s=20,c=z, marker = 'o', cmap = cm.jet )  # use a simple cmap argument!

plt.show()