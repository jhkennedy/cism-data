#!/usr/bin/env python

import os
import sys
import math
import numpy
import argparse

from netCDF4 import Dataset
from shapely.geometry import shape 

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

    parser.add_argument('-i', '--input', 
            type=abs_existing_file,
            help='NetCDF dataset to build a SCRIP grid for.')
    parser.add_argument('-d', '--dims', nargs=2,
            default=['y1','x1'],
            help='The name of the dimensions used to describe the 2D regualar '+ 
                 'grid axes, in the order they apear in the NetCDF file\'s '+
                 'varaible description.')
    parser.add_argument('-p', '--projection',
            default='epsg', type=str.lower,
            choices=['epsg','bamber'],
            help='The projection of the NetCDF dataset.')

    volume = parser.add_mutually_exclusive_group()
    volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
    volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

    return parser.parse_args(args)


def main(args):
    #FIXME: Pick projection based on input file!
    proj_epsg3413, proj_eigen_gl04c = projections.greenland()
    if args.projection == 'epsg':
        proj = proj_epsg3413
        scrip_title = "CISM EPSG:3413 Grid"
    elif args.projection == 'bamber':
        proj = proj_eigen_gl04c
        scrip_title = "CISM Bamber Grid"

    # load the dataset
    nc_base = Dataset(args.input, 'r+')
    base = projections.DataGrid()
    base.y = nc_base.variables[args.dims[0]]
    base.x = nc_base.variables[args.dims[1]]
    base.dy = base.y[1]-base.y[0]
    base.dx = base.x[1]-base.x[0]
    base.ny = base.y[:].shape[0]
    base.nx = base.x[:].shape[0]
    base.N = base.ny*base.nx
    base.make_grid()

    lon_grid, lat_grid = proj(base.x_grid.ravel(), base.y_grid.ravel(), inverse=True)
    lon_grid.shape = base.x_grid.shape 
    lat_grid.shape = base.x_grid.shape
    base.lon_grid = lon_grid
    base.lat_grid  = lat_grid


    base.lat = nc_base.createVariable('lat','f4', (args.dims[0],args.dims[1]))
    base.lat[:] = base.lat_grid[:]
    base.lat.units = 'degrees'
    
    base.lon = nc_base.createVariable('lon','f4', (args.dims[0],args.dims[1]))
    base.lon[:] = base.lon_grid[:]
    base.lon.units = 'degrees'

    nc_base.close()
    

if __name__ == '__main__':
    main(parse_args())

