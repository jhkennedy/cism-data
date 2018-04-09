#!/usr/bin/env python2

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator as RGI

from util import speak
from util import projections
from util.ncfunc import get_nc_file


def abs_existing_file(file):
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        print("Error! File does not exist: \n    "+file)
        sys.exit(1)
    return file


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)   

    parser.add_argument('-s', '--show', action='store_true', 
            help='Plot the interpolted data.')
    
    return parser.parse_args(args)


BASE = 'greenland_percent_coverage.nc'
RACMO_GRID = abs_existing_file('RACMO23_masks_ZGRN11.nc')
RACMO_DATA = abs_existing_file('racmo23_GRN_monthly.t2m.1960-2005.DJF.nc') 
PRCT_ICESHEET = abs_existing_file('surfdata_0.9x1.25_simyr1850_c110725.nc')


def main(args):
    nc_rgrid = get_nc_file(RACMO_GRID,'r')
    nc_racmo = get_nc_file(RACMO_DATA,'r')
    nc_prct = get_nc_file(PRCT_ICESHEET,'r')

    rgrid = projections.DataGrid()
    racmo = projections.DataGrid()
    prct = projections.DataGrid()
    
    
    rgrid.lat_grid = nc_rgrid.variables['lat']
    rgrid.lon_grid = nc_rgrid.variables['lon']

    racmo.lat_grid = nc_racmo.variables['LAT']
    racmo.lon_grid_original = nc_racmo.variables['LON']
    racmo.lon_grid = np.mod(racmo.lon_grid_original[:,:]+360.0, 360.0)

    prct.lat_grid = nc_prct.variables['LATIXY']
    prct.lon_grid = nc_prct.variables['LONGXY']
    prct.icesheet = nc_prct.variables['PCT_GLC_ICESHEET']

    if not np.allclose(rgrid.lat_grid[:,:], racmo.lat_grid[:,:]):
        print('The lats are different.')
        print(np.amin(rgrid.lat_grid[:,:] - racmo.lat_grid[:,:]))
        print(np.amax(rgrid.lat_grid[:,:] - racmo.lat_grid[:,:]))
        print(np.mean(rgrid.lat_grid[:,:] - racmo.lat_grid[:,:]))

    if not np.allclose(rgrid.lon_grid[:,:], racmo.lon_grid[:,:]):
        print('The lons are different.')
        print(np.amin(rgrid.lon_grid[:,:] - racmo.lon_grid[:,:]))
        print(np.amax(rgrid.lon_grid[:,:] - racmo.lon_grid[:,:]))
        print(np.mean(rgrid.lon_grid[:,:] - racmo.lon_grid[:,:]))

    prct.lat = prct.lat_grid[:,0]
    prct.lon = prct.lon_grid[0,:]

    #NOTE: bound_error and fill_value are used because there is no way to
    #      adjust for periodic nature of longitude.  Just zeroing out the %
    #      works because Greenland doesn't cross the prime meridian. This would
    #      NOT work for Antarctica. 
    rgi = RGI( points=(prct.lat, prct.lon), values=prct.icesheet[:,:], 
            method='nearest', bounds_error=False, fill_value=0.0 )

    racmo.icesheet = rgi( (racmo.lat_grid[:,:].ravel(), racmo.lon_grid[:,:].ravel()) )
    racmo.icesheet.shape = racmo.lat_grid.shape

    
    if args.show:
        fig, ax = plt.subplots(1,1, figsize=(12,10), dpi=150)
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

        plt.contourf(racmo.lon_grid_original[:,:], racmo.lat_grid[:,:], racmo.icesheet[:,:])
        plt.scatter(racmo.lon_grid_original[::50,::50], racmo.lat_grid[::50,::50], c='b')
        
        plt.show()

    
    nc_base = Dataset(BASE,'w', format='NETCDF4')
    nc_base.createDimension('lsmlat', racmo.lat_grid.shape[0])
    nc_base.createDimension('lsmlon', racmo.lat_grid.shape[1])

    base = projections.DataGrid()
    
    base.prct = nc_base.createVariable('PCT_GLC_ICESHEET', 'f4', ('lsmlat','lsmlon',) )
    base.prct[:,:] = racmo.icesheet[:,:]
    base.prct.long_name = 'percent ice sheet'
    base.prct.units = 'unitless'

   
    nc_base.close()
    nc_rgrid.close()
    nc_racmo.close()
    nc_prct.close()


if __name__ == '__main__':
    main(parse_args())

