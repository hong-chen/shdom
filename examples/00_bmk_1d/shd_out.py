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


def get_data_shdom_out_txt_old(
        fname,
        verbose=True,
        ):

    headers = []
    with open(fname, 'r') as f:
        line = f.readline().strip()
        while line[0] == '!':
            headers.append(line)
            line = f.readline().strip()

    output_type = headers[16].split('=')[-1].lower().strip()

    if output_type == 'radiance':

        if verbose:
            print('Reading SHDOM output file:\n <%s> ...' % os.path.abspath(fname))
            print('╭────────────────────────────────────────────────────────────────────────────╮')

        # retrieve Nd, Nx, Ny for radiance
        #╭────────────────────────────────────────────────────────────────────────────╮#
        line = headers[17]

        index_sx = line.find('NXO')
        line = line[index_sx:]
        index_sy = line.find('NYO')
        index_sd = line.find('NDIR')

        Nx = int(line[:index_sy].split('=')[-1].strip())
        Ny = int(line[index_sy:index_sd].split('=')[-1].strip())
        Nd = int(line[index_sd:].split('=')[-1].strip())
        #╰────────────────────────────────────────────────────────────────────────────╯#

        data_ = np.loadtxt(fname, comments='!')[:, 2:]
        data = np.zeros((Nd, Nx, Ny), dtype=np.float32)
        for i in range(Nd):
            data[i, ...] = data_[:, i].reshape((Nx, Ny))

        if verbose:
            print('%s (Nd, Nx, Ny): %s' % (output_type.title(), data_.shape))
            print('╰────────────────────────────────────────────────────────────────────────────╯')

    elif output_type == 'flux':
        if verbose:
            print('Reading SHDOM output file:\n <%s> ...' % os.path.abspath(fname))
            print('╭────────────────────────────────────────────────────────────────────────────╮')

        # retrieve Nd, Nx, Ny for radiance
        #╭────────────────────────────────────────────────────────────────────────────╮#
        if 'NXO' in headers[17]:
            line = headers[17]

            index_sx = line.find('NXO')
            line = line[index_sx:]
            index_sy = line.find('NYO')

            Nx = int(line[:index_sy].split('=')[-1].strip())
            Ny = int(line[index_sy:].split('=')[-1].strip())
            Nz = 1
        else:
            line = headers[2]
            index_sx = line.find('NX')
            line = line[index_sx:]
            index_sy = line.find('NY')
            index_sz = line.find('NZ')
            index_sp = line.find('NPTS')

            Nx = int(line[:index_sy].split('=')[-1].strip())
            Ny = int(line[index_sy:index_sz].split('=')[-1].strip())
            Nz = int(line[index_sz:index_sp].split('=')[-1].strip())
        #╰────────────────────────────────────────────────────────────────────────────╯#
        vnames = headers[-1][1:].split()
        Nd = len(vnames)
        if 'X' in vnames:
            Nd -= 1
        if 'Y' in vnames:
            Nd -= 1
        if 'Z' in vnames:
            Nd -= 1

        data_ = np.loadtxt(fname, comments='!')[:, -Nd:]
        data = np.zeros((Nd, Nx, Ny, Nz), dtype=np.float32)
        for i in range(Nd):
            data[i, ...] = data_[:, i].reshape((Nx, Ny, Nz))

        if verbose:
            print('%s (Nd, Nx, Ny, Nz): %s' % (output_type.title(), data.shape))
            print('╰────────────────────────────────────────────────────────────────────────────╯')

    else:

        if verbose:
            print('Reading SHDOM output file:\n <%s> ...' % os.path.abspath(fname))
            print('╭────────────────────────────────────────────────────────────────────────────╮')

        # retrieve Nd, Nx, Ny for radiance
        #╭────────────────────────────────────────────────────────────────────────────╮#
        if 'NXO' in headers[17]:
            line = headers[17]

            index_sx = line.find('NXO')
            line = line[index_sx:]
            index_sy = line.find('NYO')

            Nx = int(line[:index_sy].split('=')[-1].strip())
            Ny = int(line[index_sy:].split('=')[-1].strip())
            Nz = 1
        else:
            line = headers[2]
            index_sx = line.find('NX')
            line = line[index_sx:]
            index_sy = line.find('NY')
            index_sz = line.find('NZ')
            index_sp = line.find('NPTS')

            Nx = int(line[:index_sy].split('=')[-1].strip())
            Ny = int(line[index_sy:index_sz].split('=')[-1].strip())
            Nz = int(line[index_sz:index_sp].split('=')[-1].strip())
        #╰────────────────────────────────────────────────────────────────────────────╯#
        vnames = headers[-1][1:].split()
        Nd = len(vnames)
        if 'X' in vnames:
            Nd -= 1
        if 'Y' in vnames:
            Nd -= 1
        if 'Z' in vnames:
            Nd -= 1

        data_ = np.loadtxt(fname, comments='!')[:, -Nd:]
        data = np.zeros((Nd, Nx, Ny, Nz), dtype=np.float32)
        for i in range(Nd):
            data[i, ...] = data_[:, i].reshape((Nx, Ny, Nz))

        if verbose:
            print('%s (Nd, Nx, Ny, Nz): %s' % (output_type.title(), data.shape))
            print('╰────────────────────────────────────────────────────────────────────────────╯')

    return data

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


def figure_check(filename='brdf_L1f.out'):

    # fname_ori = 'ori/official/%s' % filename
    fname_ori = 'tmp-data/%s' % filename
    data_ori = get_data_shdom_out_txt(fname_ori)
    fname_new = 'tmp-data/%s' % filename
    data_new = get_data_shdom_out_txt(fname_new)

    if data_ori is None or data_new is None:
        return

    for i in range(data_ori.shape[0]):
        for j in range(data_ori.shape[1]):

            data0 = data_ori[i, j, ...]
            data1 = data_new[i, j, ...]

            if data0.shape[-1] > 1:
                data0 = np.squeeze(data0).T
            else:
                data0 = np.squeeze(data0)

            if data1.shape[-1] > 1:
                data1 = np.squeeze(data1).T
            else:
                data1 = np.squeeze(data1)

            # figure
            #╭────────────────────────────────────────────────────────────────────────────╮#
            plot = True
            if plot:
                plt.close('all')
                # plot1
                #╭──────────────────────────────────────────────────────────────╮#
                if (data0.ndim==1) and (data1.ndim==1):
                    fig = plt.figure(figsize=(8, 8))
                    fig.suptitle('%s [iData-%d, iDir-%d]' % (filename, i, j))
                    ax1 = fig.add_subplot(111, aspect='auto')

                    cs = ax1.plot(data0, color='k', lw=3.0)
                    cs = ax1.plot(data1, color='r', lw=1.5)

                    ax2 = ax1.twinx()
                    ax2.plot((data1-data0)/data0*100.0, color='green', lw=1.0)
                    ax2.axhline(0.0, color='gray', ls='--', zorder=0)
                    ax2.set_ylim((-1.0, 1.0))

                    patches_legend = [
                                      mpatches.Patch(color='black' , label='SHDOM (Ori)'), \
                                      mpatches.Patch(color='red'   , label='SHDOM (New)'), \
                                     ]
                    ax1.legend(handles=patches_legend, loc='upper right', fontsize=16)
                elif (data0.ndim==2) and (data1.ndim==2):
                    fig = plt.figure(figsize=(12, 4))
                    fig.suptitle('%s [%d]' % (filename, i))
                    ax1 = fig.add_subplot(131, aspect='auto')
                    cs = ax1.imshow(data0, origin='lower', cmap='jet', zorder=0, aspect='auto')
                    ax1.axis('off')
                    ax1.set_title('SHDOM (Ori)')

                    ax2 = fig.add_subplot(133, aspect='auto')
                    cs = ax2.imshow(data1, origin='lower', cmap='jet', zorder=0, aspect='auto')
                    ax2.axis('off')
                    ax2.set_title('SHDOM (New)', color='red')

                    ax3 = fig.add_subplot(132, aspect='auto')
                    cs = ax3.imshow((data1-data0), origin='lower', cmap='seismic', zorder=0, aspect='auto', vmin=-0.1, vmax=0.1)
                    ax3.axis('off')
                #╰──────────────────────────────────────────────────────────────╯#
                # save figure
                #╭──────────────────────────────────────────────────────────────╮#
                fig.subplots_adjust(hspace=0.35, wspace=0.35)
                _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
                fname_fig = '%s_%2.2d_%2.2d.png' % (filename, i, j)
                plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
                #╰──────────────────────────────────────────────────────────────╯#
                # plt.show()
                # sys.exit()
                plt.close(fig)
                plt.clf()
            #╰────────────────────────────────────────────────────────────────────────────╯#


if __name__ == '__main__':

    filenames = [os.path.basename(fname) for fname in sorted(glob.glob('tmp-data/*shd*out*.txt'))]

    for filename in filenames:
        print(filename)
        figure_check(filename=filename)

    pass
