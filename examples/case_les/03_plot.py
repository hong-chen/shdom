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
# import cartopy.crs as ccrs
# mpl.use('Agg')


def get_data_shdom_out_rad_txt(
        fname,
        verbose=True,
        iLine=18,
        ):

    # retrieve Nd, Nx, Ny
    #╭────────────────────────────────────────────────────────────────────────────╮#
    with open(fname, 'r') as f:
        for i in range(iLine):
            line = f.readline().strip()

    index_sx = line.find('NXO')
    line = line[index_sx:]
    index_sy = line.find('NYO')
    index_sd = line.find('NDIR')

    Nx = int(line[:index_sy].split('=')[-1].strip())
    Ny = int(line[index_sy:index_sd].split('=')[-1].strip())
    Nd = int(line[index_sd:].split('=')[-1].strip())
    #╰────────────────────────────────────────────────────────────────────────────╯#

    data_ = np.loadtxt(fname, comments='!')[:, -1]
    data = data_.reshape((Nd, Nx, Ny), order='F')

    if verbose:
        print('Reading SHDOM output file: <%s> ...' % os.path.abspath(fname))
        print('╭────────────────────────────────────────────────────────────────────────────╮')
        print('Data dimension (Nd, Nx, Ny):', data.shape)
        print('╰────────────────────────────────────────────────────────────────────────────╯')

    return data


def figure_les(fname, vza_mu=0.9):

    data = get_data_shdom_out_rad_txt(fname)

    # Ny_shift = int(-(13.0 * np.tan(np.deg2rad(30.0))//0.2 + 1))
    # Ny_shift = -38

    # Ny_shift = int(-(705.0 * np.tan(np.deg2rad(45.0))//0.2 + 1))
    # print(Ny_shift)
    Ny_shift = 0


    for i in range(data.shape[0]):

        data0 = np.roll(data[i, :, :], Ny_shift, axis=1)

        # figure
        #╭────────────────────────────────────────────────────────────────────────────╮#
        plot = True
        extent = [0, 48, 0, 48]
        if plot:
            plt.close('all')
            fig = plt.figure(figsize=(8, 8))
            # plot1
            #╭──────────────────────────────────────────────────────────────╮#
            ax1 = fig.add_subplot(111, aspect='auto')
            cs = ax1.imshow(data0.T, origin='lower', cmap='Greys_r', zorder=0, aspect='auto', vmin=0.0, vmax=0.2, extent=extent)
            # cs = ax1.imshow(data_new[i, :, :].T, origin='lower', cmap='Greys_r', zorder=0, aspect='auto', vmin=0.0, vmax=0.2) #, extent=extent, vmin=0.0, vmax=0.5)
            ax1.axis('off')
            #╰──────────────────────────────────────────────────────────────╯#
            # save figure
            #╭──────────────────────────────────────────────────────────────╮#
            fig.subplots_adjust(hspace=0.35, wspace=0.35)
            _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
            fname_fig = '%s_%2.2d.png' % (_metadata_['Function'], i)
            plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
            #╰──────────────────────────────────────────────────────────────╯#
            # plt.show()
            # sys.exit()
            plt.close(fig)
            plt.clf()
        #╰────────────────────────────────────────────────────────────────────────────╯#


if __name__ == '__main__':

    fname = 'tmp-data/shdom-out_les.txt'

    figure_les(fname)

    # Z0 = 20.0
    # vza = 30.0
    # vza_mu = np.cos(np.deg2rad(vza))
    # print(vza_mu)
    # vza_rad = np.arccos(vza_mu)
    # print(np.rad2deg(vza_rad))
    pass
