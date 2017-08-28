#!/usr/bin/python

import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

x, y, z, col = np.genfromtxt("dfs/" + sys.argv[1] + ".smooth.txt", unpack=True)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
line = ax.scatter(x, y, z, c = col, s=5, edgecolors='none', depthshade=0,linewidth=0)
ax.view_init(27, -7)
ax.tick_params(direction = "in")
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
ax.grid(b=False)
ax.w_xaxis.gridlines.set_lw(1.0)
ax.w_yaxis.gridlines.set_lw(2.0)
ax.w_zaxis.gridlines.set_lw(3.0)

#ax.w_xaxis._axinfo.update({'grid' : {'color': (0, 0, 0, 1)}})
#ax.w_yaxis._axinfo.update({'grid' : {'color': (0, 0, 0, 1)}})
#ax.w_zaxis._axinfo.update({'grid' : {'color': (0, 0, 0, 1)}})

ax.xaxis._axinfo['tick']['inward_factor'] = 0.4
ax.xaxis._axinfo['tick']['outward_factor'] = 0
ax.yaxis._axinfo['tick']['inward_factor'] = 0.4
ax.yaxis._axinfo['tick']['outward_factor'] = 0
ax.zaxis._axinfo['tick']['inward_factor'] = 0.4
ax.zaxis._axinfo['tick']['outward_factor'] = 0

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.zaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels()):
    item.set_fontsize(7)

plt.savefig("pdfs/" + sys.argv[1] + ".pdf", transparent=True)

