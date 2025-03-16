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
name_tag = '00_er3t_mca'
fdir0    = er3t.common.fdir_examples
photons = 1e9
Ncpu    = 12
rcParams['font.size'] = 24
#╰────────────────────────────────────────────────────────────────────────────╯#


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

            data_shd = np.squeeze(data[i, j, ...])

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
                cs = ax1.imshow(data_shd, origin='lower', cmap='Greys_r', zorder=0, aspect='auto', vmin=0.0, vmax=0.5, extent=extent)
                ax1.axis('off')
                #╰──────────────────────────────────────────────────────────────╯#
                # save figure
                #╭──────────────────────────────────────────────────────────────╮#
                fig.subplots_adjust(hspace=0.35, wspace=0.35)
                _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
                fname_fig = '%s_%2.2d_%2.2d.png' % (_metadata_['Function'], i, j)
                plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False, pad_inches=0.0)
                #╰──────────────────────────────────────────────────────────────╯#
                # plt.show()
                # sys.exit()
                plt.close(fig)
                plt.clf()
            #╰────────────────────────────────────────────────────────────────────────────╯#



def run_mcarats(
        wavelength=650.0,
        solver='3D',
        overwrite=True,
        plot=False
        ):

    """
    Similar to test_02 but for calculating radiance fields using LES data (nadir radiance at
    the satellite altitude of 705km)

    Additionally, Mie phase function is used instead of HG that was used in test_01 - test_04

    To run this test, we will need data/00_er3t_mca/aux/les.nc
    """

    _metadata   = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
    fdir='%s/tmp-data/%s/%s' % (fdir0, name_tag, _metadata['Function'])

    if not os.path.exists(fdir):
        os.makedirs(fdir)

    # define an atmosphere object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # levels: altitude of the layer interface in km, here, levels will be 0.0, 1.0, 2.0, ...., 20.0
    levels    = np.linspace(0.0, 20.0, 21)

    # file name of the pickle file for atmosphere
    fname_atm = '%s/atm.pk' % fdir

    # atmosphere object
    atm0      = er3t.pre.atm.atm_atmmod(levels=levels, fname=fname_atm, overwrite=overwrite)

    # data can be accessed at
    #     atm0.lev['altitude']['data']
    #     atm0.lev['pressure']['data']
    #     atm0.lev['temperature']['data']
    #     atm0.lev['h2o']['data']
    #     atm0.lev['o3']['data']
    #     atm0.lev['o2']['data']
    #     atm0.lev['co2']['data']
    #     atm0.lev['ch4']['data']
    #
    #     atm0.lay['altitude']['data']
    #     atm0.lay['pressure']['data']
    #     atm0.lay['temperature']['data']
    #     atm0.lay['h2o']['data']
    #     atm0.lay['o3']['data']
    #     atm0.lay['o2']['data']
    #     atm0.lay['co2']['data']
    #     atm0.lay['ch4']['data']
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define an absorption object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # file name of the pickle file for absorption
    fname_abs = '%s/abs.pk' % fdir

    # absorption object
    abs0 = er3t.pre.abs.abs_rep(wavelength=wavelength, fname=fname_abs, target='modis', band_name='modis_aqua_b01', atm_obj=atm0, overwrite=overwrite)

    # data can be accessed at
    #     abs0.coef['wavelength']['data']
    #     abs0.coef['abso_coef']['data']
    #     abs0.coef['slit_func']['data']
    #     abs0.coef['solar']['data']
    #     abs0.coef['weight']['data']
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define an cloud object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # file name of the netcdf file
    fname_nc  = '%s/data/00_er3t_mca/aux/les.nc' % (er3t.common.fdir_examples)

    # file name of the pickle file for cloud
    fname_les = '%s/les.pk' % fdir

    # cloud object
    cld0      = er3t.pre.cld.cld_les(fname_nc=fname_nc, fname=fname_les, coarsen=[2, 2, 25], overwrite=overwrite)

    # data can be accessed at
    #     cld0.lay['x']['data']
    #     cld0.lay['y']['data']
    #     cld0.lay['nx']['data']
    #     cld0.lay['ny']['data']
    #     cld0.lay['dx']['data']
    #     cld0.lay['dy']['data']
    #     cld0.lay['altitude']['data']
    #     cld0.lay['extinction']['data']
    #     cld0.lay['temperature']['data']
    #
    #     cld0.lev['altitude']['data']
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define mca_sca object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    pha0 = er3t.pre.pha.pha_mie_wc(wavelength=wavelength)
    sca  = er3t.rtm.mca.mca_sca(pha_obj=pha0, fname='%s/mca_sca.bin' % fdir, overwrite=overwrite)
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define mcarats 1d and 3d "atmosphere", can represent aersol, cloud, atmosphere
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # inhomogeneous 3d mcarats "atmosphere"
    atm3d0  = er3t.rtm.mca.mca_atm_3d(cld_obj=cld0, atm_obj=atm0, pha_obj=pha0, fname='%s/mca_atm_3d.bin' % fdir, overwrite=overwrite)
    # data can be accessed at
    #     atm3d0.nml['Atm_nx']['data']
    #     atm3d0.nml['Atm_ny']['data']
    #     atm3d0.nml['Atm_dx']['data']
    #     atm3d0.nml['Atm_dy']['data']
    #     atm3d0.nml['Atm_nz3']['data']
    #     atm3d0.nml['Atm_iz3l']['data']
    #     atm3d0.nml['Atm_tmpa3d']['data']
    #     atm3d0.nml['Atm_abst3d']['data']
    #     atm3d0.nml['Atm_extp3d']['data']
    #     atm3d0.nml['Atm_omgp3d']['data']
    #     atm3d0.nml['Atm_apfp3d']['data']

    # homogeneous 1d mcarats "atmosphere"
    atm1d0  = er3t.rtm.mca.mca_atm_1d(atm_obj=atm0, abs_obj=abs0)
    # data can be accessed at
    #     atm1d0.nml[ig]['Atm_zgrd0']['data']
    #     atm1d0.nml[ig]['Atm_wkd0']['data']
    #     atm1d0.nml[ig]['Atm_mtprof']['data']
    #     atm1d0.nml[ig]['Atm_tmp1d']['data']
    #     atm1d0.nml[ig]['Atm_nkd']['data']
    #     atm1d0.nml[ig]['Atm_nz']['data']
    #     atm1d0.nml[ig]['Atm_ext1d']['data']
    #     atm1d0.nml[ig]['Atm_abs1d']['data']
    #     atm1d0.nml[ig]['Atm_omg1d']['data']
    #     atm1d0.nml[ig]['Atm_apf1d']['data']


    # make them into python list, can contain more than one 1d or 3d mcarats "atmosphere"
    atm_1ds   = [atm1d0]
    atm_3ds   = [atm3d0]
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define mcarats object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # run mcarats
    mca0 = er3t.rtm.mca.mcarats_ng(
            date=datetime.datetime(2017, 8, 13),
            atm_1ds=atm_1ds,
            atm_3ds=atm_3ds,
            Ng=abs0.Ng,
            target='radiance',
            surface_albedo=0.03,
            sca=sca,
            solar_zenith_angle=30.0,
            solar_azimuth_angle=45.0,
            sensor_zenith_angle=0.0,
            sensor_azimuth_angle=0.0,
            sensor_altitude=19000.0,
            fdir='%s/%4.4d/rad_%s' % (fdir, wavelength, solver.lower()),
            Nrun=3,
            photons=photons,
            weights=abs0.coef['weight']['data'],
            solver=solver,
            Ncpu=Ncpu,
            mp_mode='py',
            overwrite=overwrite
            )

    # data can be accessed at
    #     mca0.Nrun
    #     mca0.Ng
    #     mca0.nml         (Nrun, Ng), e.g., mca0.nml[0][0], namelist for the first g of the first run
    #     mca0.fnames_inp  (Nrun, Ng), e.g., mca0.fnames_inp[0][0], input file name for the first g of the first run
    #     mca0.fnames_out  (Nrun, Ng), e.g., mca0.fnames_out[0][0], output file name for the first g of the first run
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # define mcarats output object
    #╭────────────────────────────────────────────────────────────────────────────╮#
    # read mcarats output files (binary) and save the data into h5 file
    # The mode can be specified as 'all', 'mean', 'std', if 'all' is specified, the data will have last
    # dimension of number of runs
    # e.g.,
    # out0 = mca_out_ng(fname='mca-out-rad-3d_les.h5', mca_obj=mca0, abs_obj=abs0, mode='mean', squeeze=True, verbose=True, overwrite=True)
    # out0 = mca_out_ng(fname='mca-out-rad-3d_les.h5', mca_obj=mca0, abs_obj=abs0, mode='std' , squeeze=True, verbose=True, overwrite=True)
    # out0 = mca_out_ng(fname='mca-out-rad-3d_les.h5', mca_obj=mca0, abs_obj=abs0, mode='all' , squeeze=True, verbose=True, overwrite=True)

    fname_h5 = '%s/mca-out-rad-%s_%s.h5' % (fdir, solver.lower(), _metadata['Function'])
    out0 = er3t.rtm.mca.mca_out_ng(fname=fname_h5, mca_obj=mca0, abs_obj=abs0, mode='mean', squeeze=True, verbose=True, overwrite=overwrite)

    # data can be accessed at
    #     out0.data['rad']['data']
    #╰────────────────────────────────────────────────────────────────────────────╯#


    # plot
    #╭────────────────────────────────────────────────────────────────────────────╮#
    if plot:
        fname_png = '%s-%s_%s.png' % (name_tag, _metadata['Function'], solver.lower())

        fig = plt.figure(figsize=(8, 6))
        ax1 = fig.add_subplot(111)
        cs = ax1.imshow(np.transpose(out0.data['rad']['data']), cmap='Greys_r', vmin=0.0, vmax=0.3, origin='lower')
        ax1.set_xlabel('X Index')
        ax1.set_ylabel('Y Index')
        ax1.set_title('Radiance at %.2f nm (%s Mode)' % (wavelength, solver))
        plt.savefig(fname_png, bbox_inches='tight')
        plt.close(fig)
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # References
    #╭────────────────────────────────────────────────────────────────────────────╮#
    er3t.util.print_reference()
    #╰────────────────────────────────────────────────────────────────────────────╯#


def figure_diff_vs_ipa():

    # shdom
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_ipa = 'tmp-data/shdom-out_les-ipa.txt'
    fname_3d  = 'tmp-data/shdom-out_les-ckd.txt'
    data_shd_ipa = np.squeeze(get_data_shdom_out_txt(fname_ipa))
    data_shd_3d  = np.squeeze(get_data_shdom_out_txt(fname_3d))
    data_shd_diff = data_shd_3d - data_shd_ipa
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # mcarats
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_ipa = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-ipa_run_mcarats.h5'
    fname_3d  = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-3d_run_mcarats.h5'

    f = h5py.File(fname_ipa, 'r')
    data_mca_ipa = f['mean/rad'][...].T
    f.close()

    f = h5py.File(fname_3d, 'r')
    data_mca_3d  = f['mean/rad'][...].T
    f.close()

    data_mca_diff = data_mca_3d - data_mca_ipa
    #╰────────────────────────────────────────────────────────────────────────────╯#

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
        ax1.scatter(data_shd_ipa, data_shd_diff, s=3, c='r', lw=0.0, alpha=1.0)
        ax1.scatter(data_mca_ipa, data_mca_diff, s=3, c='b', lw=0.0, alpha=1.0)
        ax1.axhline(0.0, color='gray', ls='--')
        # ax1.set_xlim((0, 1))
        # ax1.set_ylim((0, 1))
        ax1.set_xlabel('IPA')
        ax1.set_ylabel('3D - IPA')
        # ax1.set_title('Plot1')
        # ax1.xaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
        # ax1.yaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
        #╰──────────────────────────────────────────────────────────────╯#
        patches_legend = [
                         mpatches.Patch(color='red'   , label='SHDOM'), \
                         mpatches.Patch(color='blue'  , label='MCARaTS'), \
                         ]
        # ax1.legend(handles=patches_legend, bbox_to_anchor=(0., 1.01, 1., .102), loc=3, ncol=len(patches_legend), mode="expand", borderaxespad=0., frameon=False, handletextpad=0.2, fontsize=14)
        ax1.legend(handles=patches_legend, loc='upper right', fontsize=16)

        # save figure
        #╭──────────────────────────────────────────────────────────────╮#
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
        fname_fig = '%s_%s.png' % (_metadata_['Date'], _metadata_['Function'],)
        plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False)
        #╰──────────────────────────────────────────────────────────────╯#
        plt.show()
        sys.exit()
        plt.close(fig)
        plt.clf()
    #╰────────────────────────────────────────────────────────────────────────────╯#

def figure_shd_vs_mca():

    # shdom
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_ipa = 'tmp-data/shdom-out_les-ipa.txt'
    fname_3d  = 'tmp-data/shdom-out_les-ckd.txt'
    data_shd_ipa = np.squeeze(get_data_shdom_out_txt(fname_ipa))
    data_shd_3d  = np.squeeze(get_data_shdom_out_txt(fname_3d))
    data_shd_diff = data_shd_3d - data_shd_ipa
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # mcarats
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_ipa = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-ipa_run_mcarats.h5'
    fname_3d  = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-3d_run_mcarats.h5'

    f = h5py.File(fname_ipa, 'r')
    data_mca_ipa = f['mean/rad'][...].T
    f.close()

    f = h5py.File(fname_3d, 'r')
    data_mca_3d  = f['mean/rad'][...].T
    f.close()

    data_mca_diff = data_mca_3d - data_mca_ipa
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # figure
    #╭────────────────────────────────────────────────────────────────────────────╮#
    plot = True
    if plot:
        plt.close('all')
        fig = plt.figure(figsize=(8, 6))
        # fig.suptitle('Figure')
        # plot1
        #╭──────────────────────────────────────────────────────────────╮#
        ax1 = fig.add_subplot(111, aspect='equal')
        ax1.scatter(data_mca_3d, data_shd_3d, s=3, c='r', lw=0.0, alpha=1.0)
        ax1.scatter(data_mca_ipa, data_shd_ipa, s=3, c='b', lw=0.0, alpha=1.0)
        ax1.plot([0.0, 1.0], [0.0, 1.0], color='gray', ls='--')
        ax1.set_xlim((0, 0.4))
        ax1.set_ylim((0, 0.4))
        ax1.set_xlabel('MCARaTS')
        ax1.set_ylabel('SHDOM')
        # ax1.xaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
        # ax1.yaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
        #╰──────────────────────────────────────────────────────────────╯#
        patches_legend = [
                         mpatches.Patch(color='red'   , label='3D'), \
                         mpatches.Patch(color='blue'  , label='IPA'), \
                         ]
        # ax1.legend(handles=patches_legend, bbox_to_anchor=(0., 1.01, 1., .102), loc=3, ncol=len(patches_legend), mode="expand", borderaxespad=0., frameon=False, handletextpad=0.2, fontsize=14)
        ax1.legend(handles=patches_legend, loc='upper right', fontsize=16)

        # save figure
        #╭──────────────────────────────────────────────────────────────╮#
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
        fname_fig = '%s_%s.png' % (_metadata_['Date'], _metadata_['Function'],)
        plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False)
        #╰──────────────────────────────────────────────────────────────╯#
        plt.show()
        sys.exit()
        plt.close(fig)
        plt.clf()
    #╰────────────────────────────────────────────────────────────────────────────╯#


def figure_lwc_shd():

    # shdom
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_shd = '/Users/hchen/Work/mygit/shdom/data/shdom-lwc_les.txt'
    lwc_shd = np.sum(np.loadtxt(fname_shd, skiprows=5)[:, -2].reshape((240, 240, 4)), axis=-1).T
    cer_shd = np.mean(np.loadtxt(fname_shd, skiprows=5)[:, -1].reshape((240, 240, 4)), axis=-1).T
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # mcarats
    #╭────────────────────────────────────────────────────────────────────────────╮#
    fname_les = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/les.pk'
    cld0 = er3t.pre.cld.cld_les(fname=fname_les, coarsen=[2, 2, 25], overwrite=False)
    cer_mca = np.mean(cld0.lay['cer']['data'], axis=-1).T
    #╰────────────────────────────────────────────────────────────────────────────╯#

    # figure
    #╭────────────────────────────────────────────────────────────────────────────╮#
    plot = True
    if plot:
        plt.close('all')
        fig = plt.figure(figsize=(12, 6))
        # fig.suptitle('Figure')
        # plot1
        #╭──────────────────────────────────────────────────────────────╮#
        ax1 = fig.add_subplot(121)
        cs = ax1.imshow(lwc_shd, origin='lower', cmap='jet', zorder=0, vmin=0.0, vmax=1.0)

        ax2 = fig.add_subplot(122)
        cs = ax2.imshow(cer_shd-cer_mca, origin='lower', cmap='seismic', zorder=0, vmin=-0.5, vmax=0.5)

        # ax1.scatter(data_shd_x, data_shd_y, s=2, c='k', lw=0.0)
        # ax1.set_xlim((0, 1))
        # ax1.set_ylim((0, 1))
        # ax1.set_xlabel('X')
        # ax1.set_ylabel('Y')
        # ax1.set_title('Plot1')
        # ax1.xaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
        # ax1.yaxis.set_major_locator(FixedLocator(np.arange(0, 100, 5)))
        #╰──────────────────────────────────────────────────────────────╯#
        # save figure
        #╭──────────────────────────────────────────────────────────────╮#
        fig.subplots_adjust(hspace=0.35, wspace=0.35)
        _metadata_ = {'Computer': os.uname()[1], 'Script': os.path.abspath(__file__), 'Function':sys._getframe().f_code.co_name, 'Date':datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}
        fname_fig = '%s.png' % (_metadata_['Function'],)
        plt.savefig(fname_fig, bbox_inches='tight', metadata=_metadata_, transparent=False)
        #╰──────────────────────────────────────────────────────────────╯#
        plt.show()
        sys.exit()
        plt.close(fig)
        plt.clf()
    #╰────────────────────────────────────────────────────────────────────────────╯#

def figure_les_comp(fname):

    data = get_data_shdom_out_txt(fname)

    # if ('les-ckd' in fname) or ('les-crk' in fname):
    #     tag_shd = '3D'
    #     tag_mca = '3D'
    #     fname_mca = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-3d_run_mcarats.h5'
    # elif 'les-nockd' in fname:
    #     tag_shd = '3D w.o. CKD'
    #     tag_mca = '3D w. CKD'
    #     fname_mca = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-3d_run_mcarats.h5'
    # elif 'les-ipa' in fname:
    #     tag_shd = 'IPA'
    #     tag_mca = 'IPA'
    #     fname_mca = '/Users/hchen/Work/mygit/er3t/examples/tmp-data/00_er3t_mca/run_mcarats/mca-out-rad-ipa_run_mcarats.h5'
    tag_shd = '3D'
    tag_mca = '3D'
    fname_mca = '/Users/hchen/Work/mygit/er3t/projects/data/02_modis_rad-sim/pre-data.h5'

    Xs = 150
    Ys = 150
    f = h5py.File(fname_mca, 'r')
    data_mca = f['mod/rad/rad_0650'][Xs:Xs+250, Ys:Ys+250].T
    f.close()
    data_mca = data_mca[25:-25, 25:-25]

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):

            data_shd = np.squeeze(data[i, j, ...])
            data_shd = data_shd[25:-25, 25:-25]

            # figure
            #╭────────────────────────────────────────────────────────────────────────────╮#
            plot = True
            extent = [0, 48, 0, 48]
            if plot:
                plt.close('all')
                fig = plt.figure(figsize=(14, 14))
                fig.suptitle(tag_shd)
                # plot1
                #╭──────────────────────────────────────────────────────────────╮#
                ax1 = fig.add_subplot(221, aspect='equal')
                cs = ax1.imshow(data_shd, origin='lower', cmap='Greys_r', zorder=0, aspect='equal', vmin=0.0, vmax=0.5, extent=extent)
                # ax1.set_title('SHDOM (%s)' % tag_shd)
                ax1.set_title('SHDOM')

                divider = make_axes_locatable(ax1)
                cax = divider.append_axes('right', '5%', pad='3%')
                cbar = fig.colorbar(cs, cax=cax)

                ax2 = fig.add_subplot(222, aspect='equal')
                cs = ax2.imshow((data_shd-data_mca), origin='lower', cmap='seismic', zorder=0, aspect='equal', vmin=-0.1, vmax=0.1, extent=extent)
                ax2.set_title('SHDOM - MODIS')

                divider = make_axes_locatable(ax2)
                cax = divider.append_axes('right', '5%', pad='3%')
                cbar = fig.colorbar(cs, cax=cax)


                ax3 = fig.add_subplot(223, aspect='equal')
                cs = ax3.imshow(data_mca, origin='lower', cmap='Greys_r', zorder=0, aspect='equal', vmin=0.0, vmax=0.5, extent=extent)
                # ax3.set_title('MCARaTS (%s)' % tag_mca)
                ax3.set_title('MODIS')

                divider = make_axes_locatable(ax3)
                cax = divider.append_axes('right', '5%', pad='3%')
                cbar = fig.colorbar(cs, cax=cax)

                ax4 = fig.add_subplot(224, aspect='equal')

                # heatmap
                #╭──────────────────────────────────────────────────────────────╮#
                data_x = data_mca.copy()
                data_y = data_shd.copy()
                xedge = np.linspace(0.0, 0.5, 51)
                yedge = np.linspace(0.0, 0.5, 51)
                heatmap, xedge, yedge = np.histogram2d(data_x.ravel(), data_y.ravel(), bins=(xedge, yedge))
                extent = [xedge[0], xedge[-1], yedge[0], yedge[-1]]
                x = (xedge[:-1]+xedge[1:]) / 2.0
                y = (yedge[:-1]+yedge[1:]) / 2.0
                xx, yy = np.meshgrid(x, y, indexing='ij')
                levels = np.linspace(0, 200, 201)
                cmap = mpl.colormaps['jet']
                cmap.set_under('white')
                cs = ax4.contourf(xx, yy, heatmap, levels, extent=extent, extend='both', cmap=cmap, zorder=0)
                ax4.plot(extent[:2], extent[2:], color='gray', ls='--', lw=1.5)
                ax4.set_xlim(extent[:2])
                ax4.set_ylim(extent[2:])
                ax4.set_title('Heatmap')
                ax4.plot([0.0, 1.0], [0.0, 1.0], color='gray', ls='--')
                ax4.set_xlabel('MCARaTS')
                ax4.set_ylabel('SHDOM')
                divider = make_axes_locatable(ax4)
                cax = divider.append_axes('right', '5%', pad='3%')
                cbar = fig.colorbar(cs, cax=cax)
                #╰──────────────────────────────────────────────────────────────╯#


                # ax1.scatter(data_mca, data_shd, s=1)
                # ax2.set_xlim((0.0, 0.3))
                # ax2.set_ylim((0.0, 0.3))
                ax1.axis('off')
                ax2.axis('off')
                ax3.axis('off')
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

    # fname = 'tmp-data/shdom-out_les-nockd.txt'
    fname = 'tmp-data/shdom-out_modis.txt'
    # fname = 'tmp-data/shdom-out_les-ipa.txt'

    figure_les(fname)
    # figure_les_comp(fname)

    # run_mcarats(solver='IPA')

    # figure_diff_vs_ipa()
    # figure_shd_vs_mca()

    # figure_lwc_shd()
    pass
