#!/usr/bin/env python
"""Plot all spherical harmonics up to n_max degrees"""

__author__ = "Ilya Oshchepkov"
__copyright__ = "Copyright 2018, Ilya Oshchepkov"
__credits__ = ["Ilya Oshchepkov"]
__license__ = "GPL3"
__version__ = "0.1"
__maintainer__ = "Ilya Oshchepkov"
__status__ = "Prototype"


import matplotlib as mpl
import matplotlib.pylab as plt
from matplotlib import cm, colors

import numpy as np
from scipy.special import sph_harm

mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

n_max = 3
cmap = "seismic"

x = np.linspace(-np.pi, np.pi, 200)
y = np.linspace(-np.pi/2, np.pi/2, 200)
xi, yi = np.meshgrid(x, y)

phi = x.copy()
phi[x < 0] = 2 * np.pi + x[x<0]
theta = np.pi/2 - y
phii, thetai = np.meshgrid(phi, theta)


fig, axs = plt.subplots(n_max + 1, n_max + 1,
                        subplot_kw=dict(projection='mollweide'),
                        figsize=(1.5*n_max*7,n_max*7))

xticks = [-(np.pi - np.pi/3), -np.pi/3, 0.0, np.pi/3, np.pi - np.pi/3]
yticks = [-np.pi/2+np.pi/12, -np.pi/3, -np.pi/6, 0.0, np.pi/6, np.pi/3, np.pi/2 - np.pi/12]

images = []
for i in range(n_max + 1):
    for j in range(i + 1):
        sh_sp = sph_harm(j, i, phii, thetai).real * (-1)**j
        images.append(axs[i, j].pcolormesh(xi, yi , sh_sp, cmap=cmap))
        axs[i, j].grid(color='k', alpha=0.25)
        axs[i, j].set_yticks(yticks, minor=False)
        axs[i, j].tick_params(axis = 'both', labelsize=n_max*7)
        axs[i, j].set_xticks(xticks, minor=False)
        axs[i, j].set_title('n = ' + str(i) + ', k = ' + str(j), fontsize=16 + n_max*7)
        axs[i, j].set_xlabel(r'$\boldsymbol \lambda$', fontsize=16 + n_max*7)
        axs[i, j].set_ylabel(r'$\boldsymbol{\varphi}$', fontsize=16 + n_max*7)

vmin = min(image.get_array().min() for image in images)
vmax = max(image.get_array().max() for image in images)
norm = colors.Normalize(vmin=vmin, vmax=vmax)
for im in images:
    im.set_norm(norm)

for ax in axs[np.triu_indices_from(axs, k=1)]:
    ax.remove()

fig.tight_layout()

fig.savefig('sph_harmonics_n_max_' + str(n_max) + '.png', dpi=100)
