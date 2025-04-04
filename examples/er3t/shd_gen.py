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
    # fname_les = '%s/les.pk' % fdir
    fname_les = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/les.pk'
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[2, 2, 25], overwrite=False)

    cer = cld0.lay['cer']['data']

    q_factor = 2.0
    const0   = 0.75*q_factor/(1000.0*1e-6)
    lwc = cld0.lay['extinction']['data']/(const0/cer) * 1000.0
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # lwc = cld0.lay['lwc']['data']*1000.0
    # cer[cer>=25.0] = 25.0
    # cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape

    with open('shdom-lwc_les.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % lwc.shape)
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % (alt0+0.5)) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                for iz in np.arange(Nz):
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
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

    with open('data/shdom-ckd_modis-aqua-b01.txt', 'w') as f:
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

        f.write('%d\n' % atm0.lay['altitude']['data'].size)

        f.write('!\n')
        f.write('! Alt [km] | Pres [mb] | Temp [K]\n')

        alt = atm0.lay['altitude']['data'][::-1]
        temp = atm0.lay['temperature']['data'][::-1]
        pres = atm0.lay['pressure']['data'][::-1]

        indices_sort = np.argsort(abs0.coef['weight']['data'])
        kabs = abs0.coef['abso_coef']['data'][::-1, indices_sort]

        for j in range(alt.size):
            # f.write('%.3f %.1f %.1f\n' % (alt[j], pres[j], temp[j]))
            f.write('%.3f\n' % (alt[j]))

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




def gen_clear_sky_1d(
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
    levels = np.linspace(0.0, 20.0, 11)
    fname_atm = '%s/atm.pk' % fdir
    atm0 = er3t.pre.atm.atm_atmmod(levels=levels, fname=fname_atm, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # define an absorption object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_abs = '%s/abs.pk' % fdir
    abs0 = er3t.pre.abs.abs_rep(wavelength=wavelength, fname=fname_abs, target='modis', band_name='modis_aqua_b01', atm_obj=atm0, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # gen_atm_file_1d(atm0)

    gen_ckd_file(atm0, abs0)


def gen_lwc_les(
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
    Nx_c = 2
    Ny_c = 2
    fname_nc  = '%s/data/00_er3t_mca/aux/les.nc' % (er3t.common.fdir_examples)
    # fname_les = '%s/les.pk' % fdir
    fname_les = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/les.pk'
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[Nx_c, Ny_c, 25], overwrite=False)

    cer = cld0.lay['cer']['data']

    q_factor = 2.0
    const0   = 0.75*q_factor/(1000.0*1e-6)
    lwc = cld0.lay['extinction']['data']/(const0/cer) * 1000.0
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # lwc = cld0.lay['lwc']['data']*1000.0
    # cer[cer>=25.0] = 25.0
    # cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape

    # with open('shdom-lwc_les.txt', 'w') as f:
    with open('shdom-lwc_les-fast.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % lwc.shape)
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % (alt0+0.5)) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                for iz in np.arange(Nz):
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return

def gen_lwc_modis(
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
    levels = np.arange(0.0, 20.1, 0.5)
    fname_atm = '%s/atm.pk' % fdir
    atm0 = er3t.pre.atm.atm_atmmod(levels=levels, fname=fname_atm, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define an absorption object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_abs = '%s/abs.pk' % fdir
    abs0 = er3t.pre.abs.abs_rep(wavelength=wavelength, fname=fname_abs, target='modis', band_name='modis_aqua_b01', atm_obj=atm0, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define an cloud object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_cld = '/Users/hchen/Work/mygit/er3t/projects/tmp-data/02_modis_rad-sim/sim-0650.0nm/3d/cld.pk'
    cld0 = er3t.pre.cld.cld_gen_cop(
            fname=fname_cld,
            overwrite=False
            )

    cer = cld0.lay['cer']['data']
    lwc = np.zeros_like(cer)

    q_factor = 2.0
    const0   = 0.75*q_factor/(1000.0*1e-6)
    lwc[cer>0.0] = cld0.lay['extinction']['data'][cer>0.0]/(const0/cer[cer>0.0]) * 1000.0

    cer[cer==0.0] = 0.5
    cer[cer>25.0] = 25.0
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # lwc = cld0.lay['lwc']['data']*1000.0
    # cer[cer>=25.0] = 25.0
    # cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape
    Nx = 150
    Ny = 150
    Xs = 150
    Ys = 150

    with open('shdom-lwc_modis.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % (Nx, Ny, Nz))
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % (alt0+0.5)) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Xs, Xs+Nx):
            for iy in np.arange(Ys, Ys+Ny):
                for iz in np.arange(Nz):
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1-Xs), (iy+1-Ys), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return

def gen_sfc_modis(
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
    # f = h5py.File('/Users/hchen/Work/01_libera/06_mod-lrt-bmk/aux/surface_albedo_land.h5', 'r')
    f = h5py.File('/Users/hchen/Work/mygit/er3t/projects/data/02_modis_rad-sim/pre-data.h5', 'r')
    f_iso = f['mod/sfc/fiso_43_0650'][...]
    f_vol = f['mod/sfc/fvol_43_0650'][...]
    f_geo = f['mod/sfc/fgeo_43_0650'][...]
    f.close()

    Nx, Ny = f_iso.shape
    Nx = 150
    Ny = 150
    Xs = 150
    Ys = 150

    # with open('shdom-sfc_land-lsrt.txt', 'w') as f:
    with open('shdom-sfc_modis.txt', 'w') as f:
        f.write('T\n')
        # f.write('%d %d %.4f %.4f\n' % (Nx, Ny, (48.0/Nx), (48.0/Ny)))
        f.write('%d %d %.4f %.4f\n' % (Nx, Ny, 0.25, 0.25))
        for ix in np.arange(Xs, Xs+Nx):
            for iy in np.arange(Ys, Ys+Ny):
                f.write('%d %d %.2f %.6f %.6f %.6f\n' % ((ix+1-Xs), (iy+1-Ys), 285.0, f_iso[ix, iy], f_vol[ix, iy], f_geo[ix, iy]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return






def convert_prpfile():

    """
    subroutine WRITE_MIE_TABLE(MIETABFILE, POLTAB, WAVELEN1, WAVELEN2, DELTAWAVE,&
        & PARTYPE, PARDENS, RINDEX,&
        & DISTFLAG, ALPHA, GAMMA, NRETAB, SRETAB, ERETAB, LOGRE,&
        & EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)
        ! Writes the table of Mie scattering properties as a function of
        ! effective radius.  There are two types of tables output:
        ! unpolarized with Legendre (Wigner d_l_00) series of a single phase
        ! matrix element (P11) or polarized with six Wigner d-function series
        ! of the Mie phase matrix elements.
        implicit none
        integer, intent(IN) :: NRETAB, MAXRANK, NRANK(NRETAB)
        logical, intent(IN) :: POLTAB, LOGRE
        complex, intent(IN) :: RINDEX
        real, intent(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
        real, intent(IN) :: PARDENS, SRETAB, ERETAB, ALPHA, GAMMA
        real, intent(IN) :: EXTINCT(NRETAB), SSALB(NRETAB)
        real, intent(IN) :: WIGCOEF(6, 0:MAXRANK, NRETAB)
        character(LEN=1), intent(IN) :: PARTYPE, DISTFLAG
        character(LEN=256), intent(IN) :: MIETABFILE
        integer :: I, J, JP, K, L, NL, NPHASEPOL
        real REFF

        open (UNIT=3, FILE=trim(adjustl(MIETABFILE)), STATUS='REPLACE')
        if (POLTAB) then
            NPHASEPOL = 6
            write (3, '(A)') '! Polarized Mie scattering table vs. effective radius (LWC=1 g/m^3)'
        else
            NPHASEPOL = 1
            write (3, '(A)') '! Mie scattering table vs. effective radius (LWC=1 g/m^3)'
        end if
        if (DELTAWAVE < 0.0) then
            write (3, '(2(1X,F8.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
        else
            write (3, '(3(1X,F8.3),A)') WAVELEN1, WAVELEN2, DELTAWAVE, '  wavelength range and averaging step (micron)'
        end if
        write (3, '(1X,F5.3,2X,A1,A)') PARDENS, PARTYPE, '   particle density (g/cm^3) and type (Water, Ice, Aerosol)'
        write (3, '(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
        if (DISTFLAG == 'L') then
            write (3, '(F7.5,A)') ALPHA, '  lognormal log standard deviation'
        else if (DISTFLAG == 'G') then
            write (3, '(F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
        else if (DISTFLAG == 'M') then
            write (3, '(2(1X,F7.5),A)') ALPHA, GAMMA, ' modified gamma size distribution shape parameters'
        end if
        write (3, '(1X,I3,2(1X,F8.3),A)') NRETAB, SRETAB, ERETAB,&
            & '  number, starting, ending effective radius'

        do I = 1, NRETAB
            if (NRETAB <= 1) then
                REFF = SRETAB
            else if (LOGRE) then
                REFF = exp((log(ERETAB) - log(SRETAB))*FLOAT(I - 1)/(NRETAB - 1) + log(SRETAB))
            else
                REFF = (ERETAB - SRETAB)*FLOAT(I - 1)/(NRETAB - 1) + SRETAB
            end if
            write (3, '(1X,F8.4,1X,E12.5,1X,F8.6,1X,I6,A)')&
                & REFF, EXTINCT(I), SSALB(I), NRANK(I), '  Reff  Ext  Alb  Nrank'
            do J = 1, NPHASEPOL
                if (POLTAB) then
                    write (3, '(I1,1X,201(1X,F10.5))') J, (WIGCOEF(J, L, I), L=0, min(NRANK(I), 200))
                else
                    write (3, '(2X,201(1X,F10.5))') (WIGCOEF(J, L, I), L=0, min(NRANK(I), 200))
                end if
                do K = 200, NRANK(I) - 1, 200
                    write (3, '(2X,200(1X,F10.5))') (WIGCOEF(J, K + L, I), L=1, min(200, NRANK(I) - K))
                end do
            end do
        end do
        close (3)
    end subroutine WRITE_MIE_TABLE


    subroutine OUTPUT_PROPERTIES(PROPFILE, NPHASEPOL, NX, NY, NZT, &
                                 DELX, DELY, ZLEVELS, TEMP, &
                                 NPHASE, MAXCOEF, NCOEF, WIGCOEF, &
                                 EXTINCT, SSALB, IPHASE)
        implicit none
        character(LEN=256), intent(IN) :: PROPFILE
        integer, intent(IN) :: NPHASEPOL, NX, NY, NZT, NPHASE, MAXCOEF
        integer, intent(IN) :: NCOEF(0:NPHASE), IPHASE(NZT, NY, NX)
        real, intent(IN) :: DELX, DELY, ZLEVELS(NZT), TEMP(NZT)
        real, intent(IN) :: WIGCOEF(NPHASEPOL, 0:MAXCOEF, 0:NPHASE)
        real, intent(IN) :: EXTINCT(NZT, NY, NX), SSALB(NZT, NY, NX)
        integer :: NPHASEOUT, I, IX, IY, IZ, J, K, L, NL
        integer, allocatable :: PHASECONV(:)

        ! Make a phase function index conversion array to output only the ones used.
        allocate (PHASECONV(0:NPHASE))
        PHASECONV(:) = 0
        NPHASEOUT = 0
        do I = 0, NPHASE
            if (any(IPHASE(:, :, :) == I)) then
                NPHASEOUT = NPHASEOUT + 1
                PHASECONV(I) = NPHASEOUT
            end if
        end do

        ! Open the property file and output the header including the
        !   tabulated phase functions
        open (UNIT=3, FILE=trim(adjustl(PROPFILE)), STATUS='REPLACE')
        if (NPHASEPOL == 1) then
            write (3, '(A)') 'Tabulated phase function property file'
        else
            write (3, '(A)') 'Polarized tabulated phase function property file'
        end if
        write (3, '(3(1X,I5))') NX, NY, NZT
        write (3, '(2(1X,F7.4),202(1X,F8.4))') DELX, DELY, (ZLEVELS(I), I=1, NZT)
        write (3, '(1X,I4)') NPHASEOUT
        do I = 0, NPHASE
            if (PHASECONV(I) > 0) then
                NL = NCOEF(I)
                do J = 1, NPHASEPOL
                    if (NPHASEPOL == 1) then
                        write (3, '(1X,I4,5X,200(1X,F9.5))') NL, (WIGCOEF(J, L, I), L=1, min(NL, 200))
                    else
                        write (3, '(1X,I1,1X,I4,5X,201(1X,F9.5))') J, NL, (WIGCOEF(J, L, I), L=0, min(NL, 200))
                    end if
                    do K = 200, NL - 1, 200
                        write (3, '(200(1X,F9.5))') (WIGCOEF(J, K + L, I), L=1, min(NL - K, 200))
                    end do
                end do
            end if
        end do

        ! Write one line of optical properties for each grid point
        do IX = 1, NX
            do IY = 1, NY
                do IZ = 1, NZT
                    write (3, '(3(1X,I4),1X,F5.1,1X,F9.4,1X,F8.6,1X,I3)') IX, IY, IZ, &
                        TEMP(IZ), EXTINCT(IZ, IY, IX), SSALB(IZ, IY, IX), PHASECONV(IPHASE(IZ, IY, IX))
                end do
            end do
        end do
        close (3)
    end subroutine OUTPUT_PROPERTIES
 


    """

    pass

def gen_lwc_les_fast(
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
    Nx_c = 5
    Ny_c = 6
    fname_nc  = '%s/data/00_er3t_mca/aux/les.nc' % (er3t.common.fdir_examples)
    # fname_les = '%s/les.pk' % fdir
    fname_les = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/example_05_rad_les_cloud_3d/les.pk'
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[Nx_c, Ny_c, 25], overwrite=False)

    cer = cld0.lay['cer']['data']

    q_factor = 2.0
    const0   = 0.75*q_factor/(1000.0*1e-6)
    lwc = cld0.lay['extinction']['data']/(const0/cer) * 1000.0
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # lwc = cld0.lay['lwc']['data']*1000.0
    # cer[cer>=25.0] = 25.0
    # cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape

    # with open('shdom-lwc_les.txt', 'w') as f:
    with open('shdom-lwc_les-fast.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % lwc.shape)
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % (alt0+0.5)) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                for iz in np.arange(Nz):
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return

def gen_sfc_land_lsrt_fast(
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
    # f = h5py.File('/Users/hchen/Work/01_libera/06_mod-lrt-bmk/aux/surface_albedo_land.h5', 'r')
    f = h5py.File('/Users/hchen/Work/mygit/er3t/projects/data/02_modis_rad-sim/pre-data.h5', 'r')
    f_iso = er3t.util.downscale(f['mod/sfc/fiso_43_0650'][...][:400, :480], [80, 96])
    f_vol = er3t.util.downscale(f['mod/sfc/fvol_43_0650'][...][:400, :480], [80, 96])
    f_geo = er3t.util.downscale(f['mod/sfc/fgeo_43_0650'][...][:400, :480], [80, 96])
    f.close()

    Nx, Ny = f_iso.shape

    # with open('shdom-sfc_land-lsrt.txt', 'w') as f:
    with open('shdom-sfc_land-lsrt-fast.txt', 'w') as f:
        f.write('T\n')
        f.write('%d %d %.4f %.4f\n' % (Nx, Ny, 0.6, 0.5))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                f.write('%d %d %.2f %.6f %.6f %.6f\n' % ((ix+1), (iy+1), 285.0, f_iso[ix, iy], f_vol[ix, iy], f_geo[ix, iy]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return


def gen_lwc_les_zero(
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
    Nx_c = 2
    Ny_c = 2
    fname_nc  = '%s/data/00_er3t_mca/aux/les.nc' % (er3t.common.fdir_examples)
    # fname_les = '%s/les.pk' % fdir
    fname_les = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/les.pk'
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[Nx_c, Ny_c, 25], overwrite=False)

    cer = cld0.lay['cer']['data']

    q_factor = 2.0
    const0   = 0.75*q_factor/(1000.0*1e-6)
    lwc = cld0.lay['extinction']['data']/(const0/cer) * 1000.0
    lwc[...] = 0.0
    cer[...] = 0.5
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # lwc = cld0.lay['lwc']['data']*1000.0
    # cer[cer>=25.0] = 25.0
    # cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape

    # with open('shdom-lwc_les.txt', 'w') as f:
    with open('shdom-lwc_les-zero.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % lwc.shape)
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % (alt0+0.5)) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                for iz in np.arange(Nz):
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return

def gen_lwc_les_full(
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
    # fname_les = '%s/les.pk' % fdir
    fname_les = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/example_05_rad_les_cloud_3d/les.pk'
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[1, 1, 25], overwrite=False)

    cer = cld0.lay['cer']['data']

    q_factor = 2.0
    const0   = 0.75*q_factor/(1000.0*1e-6)
    lwc = cld0.lay['extinction']['data']/(const0/cer) * 1000.0
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # generate les_lwc.txt
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # lwc = cld0.lay['lwc']['data']*1000.0
    # cer[cer>=25.0] = 25.0
    # cer[cer<=0.5] = 0.5
    temp = cld0.lay['temperature']['data']

    Nx, Ny, Nz = lwc.shape

    with open('shdom-lwc_les-full.txt', 'w') as f:
        f.write('2  parameter LWC file\n')
        f.write('%d %d %d\n' % lwc.shape)
        f.write('%.4f %.4f\n' % (cld0.lay['dx']['data'], cld0.lay['dy']['data']))
        f.write('%s\n' % ' '.join([str('%.4f' % (alt0+0.5)) for alt0 in cld0.lay['altitude']['data']]))
        f.write('%s\n' % ' '.join([str('%.4f' % np.mean(temp[:, :, iz])) for iz in range(Nz)]))
        for ix in np.arange(Nx):
            for iy in np.arange(Ny):
                for iz in np.arange(Nz):
                    f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), lwc[ix, iy, iz], cer[ix, iy, iz]))
                    # f.write('%d %d %d %.6f %.2f\n' % ((ix+1), (iy+1), (iz+1), 0.0, cer[ix, iy, iz]))
    #╰────────────────────────────────────────────────────────────────────────────╯#

    return


if __name__ == '__main__':

    # gen_clear_sky_1d()
    # gen_lwc_les(overwrite=False)
    # gen_sfc_modis()
    # gen_lwc_modis()

    # gen_lwc_les_fast(overwrite=False)
    # gen_sfc_land_lsrt_fast()
    # gen_lwc_les_zero(overwrite=False)
    gen_lwc_les_full()
    pass
