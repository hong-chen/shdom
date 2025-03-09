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



def get_data_shdom_out_txt(
        fname,
        verbose=True,
        ):

    if verbose:
        print('Reading SHDOM output file:\n <%s> ...' % os.path.abspath(fname))
        print('╭────────────────────────────────────────────────────────────────────────────╮')


    headers = []
    with open(fname, 'r') as f:
        line = f.readline().strip()
        while line[0] == '!':
            headers.append(line)
            line = f.readline().strip()

    Nline_output_type = [i for i in range(len(headers)) if ('OUTPUT_TYPE=' in headers[i])][0]
    output_type = headers[Nline_output_type].split('=')[-1].lower().strip()

    # Ndata
    #╭────────────────────────────────────────────────────────────────────────────╮#
    line_xy = [line for line in headers[Nline_output_type:] if ('  X  ' in line) or ('  Y  ' in line) or (' Z ' in line)][0]
    Ndata_s = 0
    if '  X  ' in line_xy:
        Ndata_s += 1
    if '  Y  ' in line_xy:
        Ndata_s += 1
    if ' Z ' in line_xy:
        Ndata_s += 1

    data_ = np.loadtxt(fname, comments='!')[:, Ndata_s:]
    Ndata = data_.shape[-1]
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # Nx, Ny, Nz
    #╭────────────────────────────────────────────────────────────────────────────╮#
    line_xy = [line for line in headers[Nline_output_type:] if ('NXO=' in line) or ('NYO=' in line)]
    if len(line_xy) == 1:
        Nz = 1
        Nx = int(line_xy[0][1:].split('NXO=')[1].split()[0])
        Ny = int(line_xy[0][1:].split('NYO=')[1].split()[0])
        if 'NDIR' in line_xy[0]:
            Ndir = int(line_xy[0][1:].split('NDIR=')[1].split()[0])
        else:
            Ndir = 1
    elif len(line_xy) == 0:
        Nx = int(headers[2][1:].split('NX=')[1].split()[0])
        Ny = int(headers[2][1:].split('NY=')[1].split()[0])
        Nz = data_.shape[0]//Nx//Ny
        Ndir = 1
    #╰────────────────────────────────────────────────────────────────────────────╯#

    data = np.zeros((Ndata, Ndir, Nx, Ny, Nz), dtype=np.float32)
    for i in range(Ndata):
        data[i, ...] = data_[:, i].reshape((Ndir, Nx, Ny, Nz))

    if verbose:
        print('%s (Ndata, Ndir, Nx, Ny, Nz): %s' % (output_type.title(), data.shape))
        print('╰────────────────────────────────────────────────────────────────────────────╯')

    return data


def figure_les(fname):

    data = get_data_shdom_out_txt(fname)

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):

            data0 = np.squeeze(data[i, j, ...])

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
                cs = ax1.imshow(data0, origin='lower', cmap='Greys_r', zorder=0, aspect='auto', vmin=0.0, vmax=0.05, extent=extent)
                # cs = ax1.plot(data0[100, :])
                ax1.axis('off')
                #╰──────────────────────────────────────────────────────────────╯#
                # save figure
                #╭──────────────────────────────────────────────────────────────╮#
                fig.subplots_adjust(hspace=0.35, wspace=0.35)
                _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
                fname_fig = '%s_%2.2d_%2.2d.png' % (_metadata_['Function'], i, j)
                plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
                #╰──────────────────────────────────────────────────────────────╯#
                plt.show()
                sys.exit()
                plt.close(fig)
                plt.clf()
            #╰────────────────────────────────────────────────────────────────────────────╯#


if __name__ == '__main__':

    fname = 'tmp-data/shdom-out_land.txt'

    # os.system('diff data/shdom-out_les.txt %s' % fname)

    figure_les(fname)

    pass
