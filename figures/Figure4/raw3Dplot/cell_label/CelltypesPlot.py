#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

x, y, z = np.genfromtxt("PC123.txt", unpack=True)
with open('cellTypes.txt') as f:
    rgb = f.read().splitlines()

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
line = ax.scatter(x, y*-1, z, c = rgb, s=5, edgecolors='none', depthshade=0)
ax.view_init(27, -7)
ax.tick_params(direction = "in")
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')

ax.w_xaxis.gridlines.set_lw(1.0)
ax.w_yaxis.gridlines.set_lw(2.0)
ax.w_zaxis.gridlines.set_lw(3.0)

ax.xaxis._axinfo['tick']['inward_factor'] = 0.4
ax.xaxis._axinfo['tick']['outward_factor'] = 0
ax.yaxis._axinfo['tick']['inward_factor'] = 0.4
ax.yaxis._axinfo['tick']['outward_factor'] = 0
ax.zaxis._axinfo['tick']['inward_factor'] = 0.4
ax.zaxis._axinfo['tick']['outward_factor'] = 0
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))


for item in ([ax.title, ax.xaxis.label, ax.yaxis.label, ax.zaxis.label] +
	ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels()):
    item.set_fontsize(7)

plt.savefig('celltypes.pdf', transparent=True)


