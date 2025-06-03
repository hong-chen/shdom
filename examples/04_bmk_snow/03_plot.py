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


def figure_les(fname):

    data = er3t.rtm.shd.get_shd_data_out(fname)
    print(np.squeeze(data))
    print(data.shape)
    sys.exit()

    for j in range(data.shape[-1]):
        for i in range(data.shape[-2]):

            data_shd = np.squeeze(data[..., i, j])

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
                cs = ax1.imshow(data_shd.T, origin='lower', cmap='viridis', zorder=0, aspect='auto', vmin=0.0, vmax=0.15, extent=extent)
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


if __name__ == '__main__':

    fname = 'tmp-data/shdom-out_snow.txt'

    figure_les(fname)

    pass
