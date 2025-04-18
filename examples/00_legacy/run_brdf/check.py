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


import er3t

_plot_ = True

def figure_check(filename='brdf_L1f.out', plot=False):

    fname_ori = 'ori/%s' % filename
    data_ori = er3t.rtm.shd.get_shd_data_out_ori(fname_ori)
    fname_new = 'new/%s' % filename
    data_new = er3t.rtm.shd.get_shd_data_out(fname_new)

    if data_ori is None or data_new is None:
        return

    for j in range(data_ori.shape[-1]):
        for i in range(data_ori.shape[-2]):

            data0 = data_ori[..., i, j]
            data1 = data_new[..., i, j]

            if data0.shape[-1] > 1:
                data0 = np.squeeze(data0)
            else:
                data0 = np.squeeze(data0)

            if data1.shape[-1] > 1:
                data1 = np.squeeze(data1)
            else:
                data1 = np.squeeze(data1)

            print('[iVar=%d, iSet=%d] Sum(Diff) = %.4f' % (j, i, (data1-data0).sum()))

            # figure
            #╭────────────────────────────────────────────────────────────────────────────╮#
            if plot:
                plt.close('all')
                # plot1
                #╭──────────────────────────────────────────────────────────────╮#
                if (data0.ndim==1) and (data1.ndim==1):
                    fig = plt.figure(figsize=(8, 8))
                    fig.suptitle('%s [iVar-%d, iSet-%d]' % (filename, i, j))
                    ax1 = fig.add_subplot(111, aspect='auto')

                    cs = ax1.plot(data0, color='k', lw=3.0)
                    cs = ax1.plot(data1, color='r', lw=1.5)
                    patches_legend = [
                                      mpatches.Patch(color='black' , label='SHDOM (Ori)'), \
                                      mpatches.Patch(color='red'   , label='SHDOM (New)'), \
                                     ]
                    ax1.legend(handles=patches_legend, loc='upper left', fontsize=16)

                    ax2 = ax1.twinx()
                    ax2.plot((data1-data0)/data0*100.0, color='green', lw=1.0)
                    ax2.axhline(0.0, color='gray', ls='--', zorder=0)
                    ax2.set_ylim((-1.0, 1.0))

                elif (data0.ndim==2) and (data1.ndim==2):
                    fig = plt.figure(figsize=(12, 4))
                    fig.suptitle('%s [%d]' % (filename, i))
                    ax1 = fig.add_subplot(131, aspect='auto')
                    cs = ax1.imshow(data0.T, origin='lower', cmap='jet', zorder=0, aspect='auto')
                    ax1.axis('off')
                    ax1.set_title('SHDOM (Ori)')

                    ax2 = fig.add_subplot(133, aspect='auto')
                    cs = ax2.imshow(data1.T, origin='lower', cmap='jet', zorder=0, aspect='auto')
                    ax2.axis('off')
                    ax2.set_title('SHDOM (New)', color='red')

                    ax3 = fig.add_subplot(132, aspect='auto')
                    cs = ax3.imshow((data1-data0).T, origin='lower', cmap='seismic', zorder=0, aspect='auto', vmin=-0.1, vmax=0.1)
                    ax3.axis('off')
                #╰──────────────────────────────────────────────────────────────╯#
                # save figure
                #╭──────────────────────────────────────────────────────────────╮#
                fig.subplots_adjust(hspace=0.35, wspace=0.35)
                _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
                fname_fig = '%s_iset-%2.2d_ivar-%2.2d.png' % (filename, i, j)
                plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
                #╰──────────────────────────────────────────────────────────────╯#
                # plt.show()
                # sys.exit()
                plt.close(fig)
                plt.clf()
            #╰────────────────────────────────────────────────────────────────────────────╯#


if __name__ == '__main__':

    filenames = [os.path.basename(fname) for fname in sorted(glob.glob('new/*.out'))]

    for filename in filenames:
        print(filename)
        figure_check(filename=filename, plot=_plot_)

    pass
