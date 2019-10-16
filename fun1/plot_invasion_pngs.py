#
# plotting script tailored for the ECM-invasion project.
#
import os
import sys
import xml.etree.ElementTree as ET
import pathlib
import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as mplc
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from pyMCDS import pyMCDS

print('# args=',len(sys.argv))
print('args=',sys.argv)
substrate_name = 'oxygen'
colormap_name = 'viridis'
fill_circles = 1
if len(sys.argv) > 1:
    tree = ET.parse(sys.argv[1])
    xml_root = tree.getroot()
    substrate_name = xml_root.find(".//substrate").text
    colormap_name = xml_root.find(".//colormap").text  # https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    fill_circles = int(xml_root.find(".//fill_circles").text)
print("colormap = ",colormap_name)
print("substrate = ",substrate_name)
print("fill_circles = ",fill_circles)

axes_min = -750.0
axes_max = 750.   # TODO: get from input file


# TODO: make figsize a function of plot_size? What about non-square plots?
# self.fig = plt.figure(figsize=(9, 9))
# self.fig = plt.figure(figsize=(18, 18))
# self.fig = plt.figure(figsize=(15, 15))  # 
#fig = plt.figure(figsize=(9, 9))  # 
fig = plt.figure(figsize=(8.4, 7))  # 
cbar = None

def circles(x, y, s, c='b', vmin=None, vmax=None, **kwargs):
    """
    See https://gist.github.com/syrte/592a062c562cd2a98a83 

    Make a scatter plot of circles. 
    Similar to plt.scatter, but the size of circles are in data scale.
    Parameters
    ----------
    x, y : scalar or array_like, shape (n, )
        Input data
    s : scalar or array_like, shape (n, ) 
        Radius of circles.
    c : color or sequence of color, optional, default : 'b'
        `c` can be a single color format string, or a sequence of color
        specifications of length `N`, or a sequence of `N` numbers to be
        mapped to colors using the `cmap` and `norm` specified via kwargs.
        Note that `c` should not be a single numeric RGB or RGBA sequence 
        because that is indistinguishable from an array of values
        to be colormapped. (If you insist, use `color` instead.)  
        `c` can be a 2-D array in which the rows are RGB or RGBA, however. 
    vmin, vmax : scalar, optional, default: None
        `vmin` and `vmax` are used in conjunction with `norm` to normalize
        luminance data.  If either are `None`, the min and max of the
        color array is used.
    kwargs : `~matplotlib.collections.Collection` properties
        Eg. alpha, edgecolor(ec), facecolor(fc), linewidth(lw), linestyle(ls), 
        norm, cmap, transform, etc.
    Returns
    -------
    paths : `~matplotlib.collections.PathCollection`
    Examples
    --------
    a = np.arange(11)
    circles(a, a, s=a*0.2, c=a, alpha=0.5, ec='none')
    plt.colorbar()
    License
    --------
    This code is under [The BSD 3-Clause License]
    (http://opensource.org/licenses/BSD-3-Clause)
    """

    if np.isscalar(c):
        kwargs.setdefault('color', c)
        c = None

    if 'fc' in kwargs:
        kwargs.setdefault('facecolor', kwargs.pop('fc'))
    if 'ec' in kwargs:
        kwargs.setdefault('edgecolor', kwargs.pop('ec'))
    if 'ls' in kwargs:
        kwargs.setdefault('linestyle', kwargs.pop('ls'))
    if 'lw' in kwargs:
        kwargs.setdefault('linewidth', kwargs.pop('lw'))
    # You can set `facecolor` with an array for each patch,
    # while you can only set `facecolors` with a value for all.

    zipped = np.broadcast(x, y, s)
    patches = [Circle((x_, y_), s_)
            for x_, y_, s_ in zipped]
    collection = PatchCollection(patches, **kwargs)
    if c is not None:
        c = np.broadcast_to(c, zipped.shape).ravel()
        collection.set_array(c)
        collection.set_clim(vmin, vmax)

    ax = plt.gca()
    ax.add_collection(collection)
    ax.autoscale_view()
    # plt.draw_if_interactive()
    if c is not None:
        plt.sci(collection)
    # return collection

#-------------------------
def plot_cells(mcds):

    cell_ids = mcds.data['discrete_cells']['ID']
#        print(cell_ids.shape)
#        print(cell_ids[:4])

    #cell_vec = np.zeros((cell_ids.shape, 3))
    num_cells = cell_ids.shape[0]
    cell_vec = np.zeros((cell_ids.shape[0], 3))
    vec_list = ['position_x', 'position_y', 'position_z']
    for i, lab in enumerate(vec_list):
        cell_vec[:, i] = mcds.data['discrete_cells'][lab]
    xvals = cell_vec[:, 0]
    yvals = cell_vec[:, 1]
    title_str = substrate_name +', '+ str(mcds.get_time()) + " min (" + str(num_cells) + " agents)"
    #   plt.title(title_str)
    #   plt.xlim(axes_min,axes_max)
    #   plt.ylim(axes_min,axes_max)
    #   plt.scatter(xvals,yvals, s=rvals*scale_radius, c=rgbs)


    cell_vols = mcds.data['discrete_cells']['total_volume']
    cell_radii = (cell_vols* 0.75 / 3.14159)**0.3333
    #circles(xvals,yvals, s=cell_radii, fc='none')   # if you just want unfilled circles

    # In [16]: cell_df = mcds.get_cell_df()
    # In [15]: type(cell_df['cell_type'])
    # Out[15]: pandas.core.series.Series

    cell_df = mcds.get_cell_df()
    # rgb = [(x*0.5, x*0.5, abs(1.0-x)) for x in cell_df['cell_type'] ]
    rgb = [(1, 1, 1) for x in cell_df['cell_type']]
    for idx in range(len(rgb)):  # make smarter
      if cell_df['cell_type'][idx] == 1:
          rgb[idx] = (0,0,1)
    if fill_circles > 0:
        circles(xvals,yvals, s=cell_radii, fc=rgb, ec='gray')  # ec='none'  s=cell_radii     fc='red'
    else:
        circles(xvals,yvals, s=cell_radii, fc='none', ec='gray')  # ec='none'  s=cell_radii     fc='red'

    plt.xlim(axes_min, axes_max)
    plt.ylim(axes_min, axes_max)
    #   ax.grid(False)
#        axx.set_title(title_str)
    plt.title(title_str)

#-------------------------
def plot_substrate(mcds):
    global cbar
    # global current_idx, axes_max, gFileId, field_index
    # scipy.io.loadmat(full_fname, info_dict)
    # M = info_dict['multiscale_microenvironment']
    # f = M[self.field_index, :]   # 4=tumor cells field, 5=blood vessel density, 6=growth substrate
    # plt.clf()

    # oxy_conc = mcds.get_concentrations('oxygen')
    oxy_conc = mcds.get_concentrations(substrate_name)

    # Get the 2D mesh for contour plotting
    xgrid, ygrid = mcds.get_mesh(True)
 
    # We want to be able to control the number of contour levels so we
    # need to do a little set up
    # num_levels = 21
    # min_conc = plane_oxy.min()
    # max_conc = plane_oxy.max()
    # my_levels = np.linspace(min_conc, max_conc, num_levels)

    numx = numy = 75
    num_contours = 10
    vmin = 30.
    vmax = 38.

    levels = MaxNLocator(nbins=30).tick_values(vmin, vmax)
#    cmap = plt.get_cmap('PiYG')
#    cmap = plt.get_cmap('viridis')
    cmap = plt.get_cmap(colormap_name)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#    my_plot = plt.contourf(xvec,xvec,M[field_idx,:].reshape(N,N), num_contours, cmap='viridis') #'viridis'
    fix_cmap = 0
    if fix_cmap > 0:
      # my_plot = plt.contourf(xvec,xvec,M[field_idx,:].reshape(N,N), levels=levels, cmap=cmap)
      my_plot = plt.contourf(xgrid, ygrid, oxy_conc.reshape(numy, numx), levels=levels, extend='both', cmap=cmap)
    else:
      # my_plot = plt.contourf(xvec,xvec,M[field_idx,:].reshape(N,N), cmap=cmap)
      my_plot = plt.contourf(xgrid, ygrid, oxy_conc.reshape(numy, numx), cmap=cmap)

    if cbar == None:
#      cbar = plt.colorbar(my_plot, boundaries=np.arange(vmin, vmax, 1.0))
      cbar = plt.colorbar(my_plot)
    else:
      cbar = plt.colorbar(my_plot, cax=cbar.ax)

    plt.xlim(-750,750)  # axes_min, axes_max)
    plt.ylim(-750,750)  # axes_min, axes_max)

#for frame in range(1180,1185):
for frame in range(0,800,100):
    fname = "output%08d" % frame
    fname_xml = "output%08d.xml" % frame
    print(fname_xml)
# full_fname = os.path.join(self.output_dir, fname)
#    full_fname = fname
#    mcds = pyMCDS(fname, self.output_dir)
    mcds = pyMCDS(fname_xml)
    print(mcds.get_substrate_names())
    # print(mcds.get_time())
    plot_substrate(mcds)
    plot_cells(mcds)
    #plt.show()
    png_file = "plot%08d.png" % frame
    fig.savefig(png_file)

