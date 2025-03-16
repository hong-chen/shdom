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


def cal_sol_fac(dtime):

    """
    er3t.util.cal_sol_fac
    Calculate solar factor that accounts for Sun-Earth distance
    Input:
        dtime: datetime.datetime object
    Output:
        solfac: solar factor
    """

    doy = dtime.timetuple().tm_yday
    eps = 0.0167086
    perh= 4.0
    rsun = (1.0 - eps*np.cos(0.017202124161707175*(doy-perh)))
    solfac = 1.0/(rsun**2)

    return solfac


def cal_shd_saa(normal_azimuth_angle):

    """
    Convert normal azimuth angle (0 pointing north, positive when clockwise) to viewing azimuth in SHDOM

    Input:
        normal_azimuth_angle: float/integer, normal azimuth angle (0 pointing north, positive when clockwise)

    Output:
        SHDOM solar azimuth angle (0 sun shining from west, positive when counterclockwise)
    """

    while normal_azimuth_angle < 0.0:
        normal_azimuth_angle += 360.0

    while normal_azimuth_angle > 360.0:
        normal_azimuth_angle -= 360.0

    shd_azimuth = 270.0 - normal_azimuth_angle
    if shd_azimuth < 0.0:
        shd_azimuth += 360.0

    return shd_azimuth


def cal_shd_vaa(normal_azimuth_angle):

    """
    Convert normal azimuth angle (0 pointing north, positive when clockwise) to viewing azimuth in SHDOM

    Input:
        normal_azimuth_angle: float/integer, normal azimuth angle (0 pointing north, positive when clockwise)

    Output:
        SHDOM sensor azimuth angle (0 sensor looking from east, positive when counterclockwise)
    """

    while normal_azimuth_angle < 0.0:
        normal_azimuth_angle += 360.0

    while normal_azimuth_angle > 360.0:
        normal_azimuth_angle -= 360.0

    shd_azimuth = 90.0 - normal_azimuth_angle
    if shd_azimuth < 0.0:
        shd_azimuth += 360.0

    return shd_azimuth


if __name__ == '__main__':

    date = datetime.datetime(2019, 9, 2)
    sza = 34.93346840064262
    saa = -144.90807471473198
    vza = 14.44008093613803
    vaa = -99.99851773953723


    print('SOLARFLUX: %.6f' % cal_sol_fac(date))
    print('SOLARMU: %.6f' % np.cos(np.deg2rad(sza)))
    print('SOLARAZ: %.6f' % cal_shd_saa(saa))

    print('SENSORMU: %.6f' % np.cos(np.deg2rad(vza)))
    print('SENSORAZ: %.6f' % cal_shd_vaa(vaa))
    pass
