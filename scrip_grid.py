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

"""
Create a SCRIP formatted grid. 
"""

SURFACE_AREA_EARTH = 510065621724000.0 # unit: m^2
#source: http://home.vikenfiber.no/humrum/WGS84_Eng.html
SQR_DEG_ON_SPHERE = 41252.96   # unit: deg^2 
SQR_RAD_ON_SPHERE = 4.*math.pi # unit rad^2

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
    nc_base = Dataset(args.input, 'r')
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

    base.ll_y = base.y_grid.flatten(order='C') - base.dy/2.  
    base.ll_x = base.x_grid.flatten(order='C') - base.dx/2. 
    base.lr_y = base.y_grid.flatten(order='C') - base.dy/2.  
    base.lr_x = base.x_grid.flatten(order='C') + base.dx/2. 
    base.ur_y = base.y_grid.flatten(order='C') + base.dy/2.  
    base.ur_x = base.x_grid.flatten(order='C') + base.dx/2. 
    base.ul_y = base.y_grid.flatten(order='C') + base.dy/2.  
    base.ul_x = base.x_grid.flatten(order='C') - base.dx/2. 

    base.ll_lon, base.ll_lat = proj(base.ll_x, base.ll_y, inverse=True)
    base.lr_lon, base.lr_lat = proj(base.lr_x, base.lr_y, inverse=True)
    base.ur_lon, base.ur_lat = proj(base.ur_x, base.ur_y, inverse=True)
    base.ul_lon, base.ul_lat = proj(base.ul_x, base.ul_y, inverse=True)

    base.corner_lat = numpy.column_stack((base.ll_lat, base.lr_lat, base.ur_lat, base.ul_lat))
    base.corner_lon = numpy.column_stack((base.ll_lon, base.lr_lon, base.ur_lon, base.ul_lon))

    min_lat = numpy.amin(base.corner_lat)
    max_lat = numpy.amax(base.corner_lat)

    min_lon = numpy.amin(base.corner_lon)
    max_lon = numpy.amax(base.corner_lon)

    proj_aea = projections.equal_area(min_lat, max_lat, (max_lon+min_lon)/2.)

    # get the area for each grid cell
    sys.stdout.write("   [%-60s] %d%%" % ('='*0, 0.))
    sys.stdout.flush()
    base.area = numpy.zeros(base.N)
    for ii in range(base.N):
        ctr = (ii*60)/base.N
        if not (ii % 100): 
            sys.stdout.write("\r   [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
            sys.stdout.flush()
        
        lat = base.corner_lat[ii,:]
        lon = base.corner_lon[ii,:]
        x, y = proj_aea(lon, lat)
       
        points = {'type':'polygon', 'coordinates':[zip(x,y)]}
        base.area[ii] = shape(points).area #m^2

    sys.stdout.write("\r   [%-60s] %d%%\n" % ('='*60, 100.))

    path_scrip, name_scrip = os.path.split(args.input)
    lc_scrip = os.path.join(path_scrip, 'SCRIPgrid_'+name_scrip)

    nc_scrip = Dataset(lc_scrip, 'w', format='NETCDF4') 
    nc_scrip.createDimension('grid_size', base.N)
    nc_scrip.createDimension('grid_corners', 4)
    nc_scrip.createDimension('grid_rank', 2)
    nc_scrip.title = scrip_title 
    nc_scrip.note = 'Created by Joseph H. Kennedy using pyproj.'

    scrip = projections.DataGrid()
    scrip.dims = nc_scrip.createVariable('grid_dims','i4', ('grid_rank'))
    #NOTE: SCRIP is a Fortran program and as such assumes data is stored in column-major order. 
    #      Because (1) netCDF stores data 'C-style' with row-major order and (2) the scrip grid 
    #      formate uses flattened arrays, it's easiest to flatten the data 'C-syle', which is the 
    #      default for numpy, and flip the dimensions. SCRIP will then read in the array structure
    #      correctly, and the expected ordering within the netCDF file we are creating will be 
    #      maintained.
    #scrip.dims[:] = numpy.array(base.dims)[::-1]
    #NOTE: Now dumping everything 'F-style' (column-major) order. 
    scrip.dims[:] = numpy.array(base.dims)[::]

    scrip.imask = nc_scrip.createVariable('grid_imask','i4', ('grid_size'))
    scrip.imask[:] = numpy.ones(base.N) 
    scrip.imask.units = 'unitless'

    scrip.center_lat = nc_scrip.createVariable('grid_center_lat','f4', ('grid_size'))
    scrip.center_lat[:] = base.lat_grid[:,:].flatten(order='C')
    scrip.center_lat.setncattr('units', 'degrees')

    scrip.center_lon = nc_scrip.createVariable('grid_center_lon','f4', ('grid_size'))
    scrip.center_lon[:] = base.lon_grid[:,:].flatten(order='C')
    scrip.center_lon.setncattr('units', 'degrees')

    scrip.corner_lat = nc_scrip.createVariable('grid_corner_lat','f4', ('grid_size','grid_corners',))
    scrip.corner_lat[:,:] = base.corner_lat[:,:]
    scrip.corner_lat.units = 'degrees'

    scrip.corner_lon = nc_scrip.createVariable('grid_corner_lon','f4', ('grid_size','grid_corners',))
    scrip.corner_lon[:,:] = base.corner_lon[:,:]
    scrip.corner_lon.units = 'degrees'

    scrip.area = nc_scrip.createVariable('grid_area','f4', ('grid_size',))
    scrip.area[:] = (base.area[:]/SURFACE_AREA_EARTH)*SQR_DEG_ON_SPHERE
    scrip.area.units = 'square degrees'

    nc_base.close()
    nc_scrip.close()
    os.chmod(lc_scrip, 0o644)   # uses an octal number!

if __name__ == '__main__':
    main(parse_args())

