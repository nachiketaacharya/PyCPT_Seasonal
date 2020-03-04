import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Make some illustrative fake data:

x = np.arange(0, np.pi, 0.1)
y = np.arange(0, 2 * np.pi, 0.1)
X, Y = np.meshgrid(x, y)
Z = np.cos(X) * np.sin(Y) * 10

colors = [(238, 43, 51),
(255, 57, 67),
(253, 123, 91),
(248, 175, 123),
(254, 214, 158),
(252, 239, 188),
(255, 254, 241),
(244, 255,255),
(187, 252, 255),
(160, 235, 255),
(123, 210, 255),
(89, 179, 238),
(63, 136, 254),
(52, 86, 254)
]
colors = [ (colors[i][0] / 255.0, colors[i][1] / 255.0, colors[i][2] / 255.0) for i in range(len(colors))]
colors.reverse()
n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
cmap_name = 'my_list'
fig, ax = plt.subplots(1, 1, figsize=(6, 9))

# Create the colormap
cm = LinearSegmentedColormap.from_list(
    cmap_name, colors, N=len(colors))
im = ax.pcolormesh(X, Y, Z, cmap=cm)

# Fewer bins will result in "coarser" colomap interpolation
ax.set_title("N bins: %s" % len(colors))
fig.colorbar(im, ax=ax)
plt.show()
