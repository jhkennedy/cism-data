#!/usr/bin/env python

import os
import sys
import scipy
import argparse

from netCDF4 import Dataset

from util import projections

"""
Create (lat,lon) grids from the (y0,x0) grids. 
"""

def abs_existing_file(file):
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        print("Error! File does not exist: \n    "+file)
        sys.exit(1)
    return file

# parse the command line arguments
parser = argparse.ArgumentParser()   # -h or --help automatically included!

parser.add_argument('-i', '--input', 
        type=abs_existing_file,
        help='netCDF dataset to build a (lat0,lon0) grid from the (y0,x0) grid.')

parser.add_argument('-p', '--projection',
        default='epsg', type=str.lower,
        choices=['epsg','bamber'],
        help='The projection of the NetCDF dataset.')

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

# get the projections
proj_epsg3413, proj_eigen_gl04c = projections.greenland()
#FIXME: Pick projection based on input file!
if args.projection == 'epsg':
    proj = proj_epsg3413
    scrip_title = "CISM EPSG:3413 Grid"
elif args.projection == 'bamber':
    proj = proj_eigen_gl04c
    scrip_title = "CISM Bamber Grid"

# load the dataset
nc_base = Dataset(args.input, 'r')
base = projections.DataGrid()
base.y = nc_base.variables['y1']
base.x = nc_base.variables['x1']
base.dy = base.y[1]-base.y[0]
base.dx = base.x[1]-base.x[0]
base.ny = base.y[:].shape[0]
base.nx = base.x[:].shape[0]
base.make_grid()

base.y0 = (base.y[:-1] + base.y[1:])/2.
base.x0 = (base.x[:-1] + base.x[1:])/2.
base.y0_grid, base.x0_grid = scipy.meshgrid(base.y0[:], base.x0[:], indexing='ij')
base.dim0 = (base.ny-1, base.nx-1)

lon0_grid, lat0_grid = proj(base.x0_grid.ravel(), base.y0_grid.ravel(), inverse=True)
lon0_grid.shape = base.x0_grid.shape 
lat0_grid.shape = base.x0_grid.shape

base.lon0_grid = lon0_grid
base.lat0_grid = lat0_grid


path_grid0, name_grid0 = os.path.split(args.input)
lc_grid0 = os.path.join(path_grid0, 'grid0_'+name_grid0)

nc_grid0 = Dataset(lc_grid0, 'w', format='NETCDF4')
nc_grid0.createDimension('time', None )
nc_grid0.createDimension('y0', base.ny-1)
nc_grid0.createDimension('x0', base.nx-1)

grid0 = projections.DataGrid()

grid0.time = nc_grid0.createVariable('time', 'f4', ('time',))
grid0.time.long_name = "Model time"
grid0.time.standard_name = "time"
grid0.time.units = "common_year since 1-1-1 0:0:0"
grid0.time.calendar = "none"
grid0.time[0] = 0.

grid0.y0   = nc_grid0.createVariable('y0',   'f4', ('y0',)  )
grid0.y0.long_name = "Cartesian y-coordinate, velocity grid"
grid0.y0.units = "meter"
grid0.y0.axis = "Y"
grid0.y0[:] = base.y0[:]

grid0.x0   = nc_grid0.createVariable('x0',   'f4', ('x0',)  )
grid0.x0.long_name = "Cartesian x-coordinate, velocity grid"
grid0.x0.units = "meter"
grid0.x0.axis = "X"
grid0.x0[:] = base.x0[:]

grid0.lon0_grid = nc_grid0.createVariable('lon0', 'f4', ('time','y0','x0',))
grid0.lon0_grid.units = 'degrees'
grid0.lon0_grid.long_name = 'grid center longitude, velocity grid'
grid0.lon0_grid.source = 'Joseph H. Kennedy, ORNL'
grid0.lon0_grid[0,:,:] = base.lon0_grid[:,:]

grid0.lat0_grid = nc_grid0.createVariable('lat0', 'f4', ('time','y0','x0',))
grid0.lat0_grid.units = 'degrees'
grid0.lat0_grid.lat_name = 'grid center latitude, velocity grid'
grid0.lat0_grid.source = 'Joseph H. Kennedy, ORNL'
grid0.lat0_grid[0,:,:] = base.lat0_grid[:,:]

nc_base.close()
nc_grid0.close()
os.chmod(lc_grid0, 0o644)   # uses an octal number!

