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

def read_file(fname):

    data = np.loadtxt(fname, comments='!')[:, -1]

    # data = data.reshape((352, -1, 352))
    # data = data.reshape((352, 352, -1))
    # data = data.reshape((-1, 352, 352))
    data = data.reshape((-1, 240, 240))
    # data = data.reshape((-1, 480, 480))

    for i in range(data.shape[0]):
        # figure
        #╭────────────────────────────────────────────────────────────────────────────╮#
        plot = True
        if plot:
            plt.close('all')
            fig = plt.figure(figsize=(8, 6))
            # fig.suptitle('Figure')
            # plot1
            #╭──────────────────────────────────────────────────────────────╮#
            ax1 = fig.add_subplot(111)
            # cs = ax1.imshow(data[i, 20:-20, 20:-20], origin='lower', cmap='Greys_r', zorder=0) #, extent=extent, vmin=0.0, vmax=0.5)
            cs = ax1.imshow(data[i, :, :], origin='lower', cmap='Greys_r', zorder=0) #, extent=extent, vmin=0.0, vmax=0.5)
            # ax1.scatter(data_x, data_y, s=2, c='k', lw=0.0)
            # ax1.set_xlim((0, 1))
            # ax1.set_ylim((0, 1))
            # ax1.set_xlabel('X')
            # ax1.set_ylabel('Y')
            # ax1.set_title('Plot1')
            # ax1.xaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
            # ax1.yaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
            ax1.axis('off')
            #╰──────────────────────────────────────────────────────────────╯#
            # save figure
            #╭──────────────────────────────────────────────────────────────╮#
            fig.subplots_adjust(hspace=0.35, wspace=0.35)
            _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
            fname_fig = '%s_%2.2d.png' % (_metadata_['Function'], i)
            plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False)
            #╰──────────────────────────────────────────────────────────────╯#
            # plt.show()
            # sys.exit()
            plt.close(fig)
            plt.clf()
        #╰────────────────────────────────────────────────────────────────────────────╯#



if __name__ == '__main__':

    fname = 'tmp-data/shdom-out_les.txt'
    read_file(fname)
    pass
