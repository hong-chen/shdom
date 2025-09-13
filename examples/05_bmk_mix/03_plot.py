import os
import sys
import glob
import datetime
import copy
import multiprocessing as mp
from collections import OrderedDict
# from tqdm import tqdm
import h5py
from pyhdf.SD import SD, SDC
from netCDF4 import Dataset
import numpy as np
from scipy import interpolate
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as mpl_path
import matplotlib.image as mpl_img
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib import rcParams, ticker
from matplotlib.ticker import FixedLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
# mpl.use('Agg')

import er3t


def figure_les(fname):

    data = er3t.rtm.shd.get_shd_data_out(fname)

    for j in range(data.shape[-1]):
        for i in range(data.shape[-2]):

            data_shd = np.squeeze(data[..., i, j])

            # figure
            #╭────────────────────────────────────────────────────────────────────────────╮#
            plot = True
            extent = [-180, 180, -90, 90]
            if plot:
                plt.close('all')
                fig = plt.figure(figsize=(12, 6))
                # plot1
                #╭──────────────────────────────────────────────────────────────╮#
                ax1 = fig.add_subplot(111)
                cs = ax1.imshow(data_shd.T, origin='lower', cmap='Greys_r', zorder=0, vmin=0.0, vmax=0.4, extent=extent)
                ax1.axis('off')
                #╰──────────────────────────────────────────────────────────────╯#
                # save figure
                #╭──────────────────────────────────────────────────────────────╮#
                fig.subplots_adjust(hspace=0.35, wspace=0.35)
                _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
                fname_fig = '%s_iset-%2.2d_ivar-%2.2d.png' % (_metadata_['Function'], i, j)
                plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
                #╰──────────────────────────────────────────────────────────────╯#
                # plt.show()
                # sys.exit()
                plt.close(fig)
                plt.clf()
            #╰────────────────────────────────────────────────────────────────────────────╯#

def figure_merra(fname):

    data = er3t.rtm.shd.get_shd_data_out(fname)
    lon_1d = np.linspace(-180.0, 180.0, 576)
    lat_1d = np.linspace(-90.0, 90.0, 361)
    lon_2d, lat_2d = np.meshgrid(lon_1d, lat_1d, indexing='ij')

    icount = 0
    for j in range(data.shape[-1]):
        for i in range(data.shape[-2]):

            data_shd = np.squeeze(data[..., i, j])

            proj = ccrs.NearsidePerspective(central_longitude=-105.24584926663005, central_latitude=40.0118810321591)
            proj0 = ccrs.PlateCarree()

            # figure
            #╭────────────────────────────────────────────────────────────────────────────╮#
            plot = True
            # extent = [-180, 180, -90, 90]
            if plot:
                plt.close('all')
                fig = plt.figure(figsize=(12, 12))
                # plot1
                #╭──────────────────────────────────────────────────────────────╮#
                ax1 = fig.add_subplot(111, projection=proj)
                ax1.pcolormesh(lon_2d, lat_2d, data_shd, vmin=0.0, vmax=0.55, cmap='Greys_r', transform=proj0)
                ax1.set_global()
                ax1.coastlines(color='gray', lw=0.5)
                g1 = ax1.gridlines()
                g1.xlocator = FixedLocator(np.arange(-180, 181, 60))
                g1.ylocator = FixedLocator(np.arange(-80, 81, 20))
                ax1.axis('off')
                #╰──────────────────────────────────────────────────────────────╯#
                # save figure
                #╭──────────────────────────────────────────────────────────────╮#
                fig.subplots_adjust(hspace=0.35, wspace=0.35)
                _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
                fname_fig = '%s_iset-%2.2d_ivar-%2.2d.png' % (_metadata_['Function'], i, j)
                plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
                #╰──────────────────────────────────────────────────────────────╯#
                # plt.show()
                # sys.exit()
                plt.close(fig)
                plt.clf()
            #╰────────────────────────────────────────────────────────────────────────────╯#
            icount += 1


if __name__ == '__main__':

    fname = 'tmp-data/shdom-out_globe.txt'
    figure_les(fname)

    # fname = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_shd/example_07_at3d_rad_cloud_merra/0550/rad_3d/shdom-out_g-000.txt'
    # figure_merra(fname)

    pass
