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

# global variables
#╭────────────────────────────────────────────────────────────────────────────╮#
name_tag = '00_er3t_shd'
fdir0    = er3t.common.fdir_examples
Ncpu    = 12
rcParams['font.size'] = 14
#╰────────────────────────────────────────────────────────────────────────────╯#

_FDIR_BIN_ = '/Users/hchen/Work/mygit/shdom/bin'
_PUT_ = '%s/put' % _FDIR_BIN_
_PROPGEN_ = '%s/propgen' % _FDIR_BIN_



def gen_two_params_lwc_file(
        wavelength=650.0,
        solver='3D',
        overwrite=False,
        plot=True
        ):

    _metadata   = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    fdir='%s/tmp-data/%s/%s' % (fdir0, name_tag, _metadata['Function'])

    if not os.path.exists(fdir):
        os.makedirs(fdir)

    # define an atmosphere object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    levels = np.linspace(0.0, 20.0, 21)
    fname_atm = '%s/atm.pk' % fdir
    atm0 = er3t.pre.atm.atm_atmmod(levels=levels, fname=fname_atm, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define an absorption object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_abs = '%s/abs.pk' % fdir
    abs0 = er3t.pre.abs.abs_16g(wavelength=wavelength, fname=fname_abs, atm_obj=atm0, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define an cloud object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_nc  = '%s/data/00_er3t_mca/aux/les.nc' % (er3t.common.fdir_examples)
    fname_les = '%s/les.pk' % fdir
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[1, 1, 25], overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    lwc = cld0.lay['lwc']['data']*1000.0
    cer = cld0.lay['cer']['data']
    cer[cer>=25.0] = 25.0
    cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape

    with open('les_lwc.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % lwc.shape)
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % alt0) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                for iz in np.arange(Nz):
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return

def gen_sfc_file(
        wavelength=650.0,
        solver='3D',
        overwrite=False,
        plot=True
        ):

    _metadata   = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    fdir='%s/tmp-data/%s/%s' % (fdir0, name_tag, _metadata['Function'])

    if not os.path.exists(fdir):
        os.makedirs(fdir)

    # generate sfc_lanb.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    f = h5py.File('/Users/hchen/Work/01_libera/06_mod-lrt-bmk/aux/surface_albedo_land.h5', 'r')
    sfc_alb = f['alb_659'][...]
    f.close()

    sfc_alb[sfc_alb<0.0] = 0.0
    sfc_alb[sfc_alb>1.0] = 1.0
    Nx, Ny = sfc_alb.shape

    with open('sfc_land.txt', 'w') as f:
        f.write('L\n')
        f.write('%d %d %.4f %.4f\n' % (Nx, Ny, (48.0/Nx), (48.0/Ny)))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                f.write('%d %d %.2f %.6f\n' % ((ix+1), (iy+1), 285.0, sfc_alb[ix, iy]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return

def gen_prop_1d(
        name_exe='propgen',
        fdir_bin=_FDIR_BIN_,
        ):

    fname_exe = os.path.join(fdir_bin, name_exe)

    command = '%s '
    pass


def gen_ckd_file(atm0, abs0):

    Nband = 1

    with open('data/modis_aqua_b01_shdom-ckd.txt', 'w') as f:
        f.write('! Correlated k-dist file: Aqua MODIS Band-1\n')
        f.write('%d ! number of bands\n' % 1)
        f.write('! Band# | Wave#1 | Wave#2 | SolFlx | Ng | g1 | g2 | ...\n')

        # wave number cm^-1
        wvln_min = 1.0/abs0.wvl_max_*1e7
        wvln_max = 1.0/abs0.wvl_min_*1e7

        sol = (abs0.coef['solar']['data']*abs0.coef['weight']['data']).sum()

        Ng = abs0.coef['weight']['data'].size

        g = ' '.join(['%.6f' % value for value in abs0.coef['weight']['data']])

        for iband in range(Nband):
            f.write('%d %.1f %.1f %.4f %d %s\n' % (iband+1, wvln_min, wvln_max, sol, Ng, g))

        f.write('%d ! Number of levels\n' % abs0.atm0.lay['altitude']['data'].size)

        f.write('!\n')
        f.write('! Alt [km] | Pres [mb] | Temp [K]\n')

        alt = atm0.lay['altitude']['data'][::-1]
        temp = atm0.lay['temperature']['data'][::-1]
        pres = atm0.lay['pressure']['data'][::-1]

        indices_sort = np.argsort(abs0.coef['weight']['data'])
        kabs = abs0.coef['abso_coef']['data'][::-1, indices_sort] * 1000.0

        for j in range(alt.size):
            f.write('%.3f %.1f %.1f\n' % (alt[j], pres[j], temp[j]))

        f.write('! iBand | iLay | AbsCoef [km^-1]\n')

        for iband in range(Nband):
            for j in range(alt.size):
                kabs_s = ' '.join(['%.5e' % kabs0 for kabs0 in kabs[j, :]])
                f.write('%d %d %s\n' % (iband+1, j+1, kabs_s))

def gen_atm_file_1d(atm0):

    fname_prp = 'data/shdom-prp_atm-clear.txt'

    parfile = 'NONE'

    maxnewphase = 2

    asymtol = 0.01

    fracphasetol = 0.1

    wavelen = 0.65

    Nzother = atm0.lay['altitude']['data'].size

    ZTother = zip(atm0.lay['altitude']['data'], atm0.lay['temperature']['data'])

    ZTother_list = []
    for ZTother0 in ZTother:
        ZTother0_s = ' '.join(['%.2f' % item for item in list(ZTother0)])
        ZTother_list.append(ZTother0_s)

    ZTother_s = '  '.join(ZTother_list)

    sfcpres = 1013.25
    #   set ZTother=(0 288  3 269  6 249  9 230  12 217  15 217  18 217 \
    #       21 218  24 221  27 224  30 227)
    polout='NONE'

    command = '%s %d %s %d %.2f %.2f %.4f \
%.2f %d %s %s %s | %s' % (_PUT_, 0, parfile, maxnewphase, asymtol, fracphasetol, wavelen,\
            sfcpres, Nzother, ZTother_s, polout, fname_prp, _PROPGEN_\
            )

    os.system(command)

    '''
    /Users/hchen/Work/mygit/shdom/bin/put 0 NONE 2 0.01 0.1  0.85 1013.25 11 0 288 3 269 6 249 9 230 12 217 15 217 18 217 21 218 24 221 27 224 30 227 P atmos_ray0.85.prp  | /Users/hchen/Work/mygit/shdom/bin/propgen| tail -3
    '''

    #   ${put_exe} 0  $parfile  $maxnewphase $asymtol $fracphasetol \
    #       $wavelen $sfcpres $Nzother $ZTother $polout $prpfile \
    #     | ${propgen_exe}| tail -3


def gen_clear_sky(
        wavelength=650.0,
        solver='3D',
        overwrite=True,
        plot=True
        ):

    _metadata   = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    fdir='%s/tmp-data/%s/%s' % (fdir0, name_tag, _metadata['Function'])

    if not os.path.exists(fdir):
        os.makedirs(fdir)

    # define an atmosphere object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    levels = np.linspace(0.0, 20.0, 21)
    fname_atm = '%s/atm.pk' % fdir
    atm0 = er3t.pre.atm.atm_atmmod(levels=levels, fname=fname_atm, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # define an absorption object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_abs = '%s/abs.pk' % fdir
    abs0 = er3t.pre.abs.abs_rep(wavelength=wavelength, fname=fname_abs, target='modis', band_name='modis_aqua_b01', atm_obj=atm0, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#

    gen_atm_file_1d(atm0)

    # gen_ckd_file(atm0, abs0)

if __name__ == '__main__':

    # gen_two_params_lwc_file(overwrite=True)
    gen_clear_sky()

    pass
