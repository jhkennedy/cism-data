#!/usr/bin/env python

import os
import sys
import math
import numpy
import argparse

from netCDF4 import Dataset
from shapely.geometry import shape 

from util import projections

"""
Create a SCRIP formatted grid. 
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
        help='netCDF dataset to build a SCRIP grid for.')

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

# get the projections
lc_bamber = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
proj_epsg3413, proj_eigen_gl04c = projections.greenland(args, lc_bamber)

surface_area_earth = 510065621724000.0 # unit: m^2
#source: http://home.vikenfiber.no/humrum/WGS84_Eng.html
sqr_deg_on_sphere = 41252.96   # unit: deg^2 
sqr_rad_on_sphere = 4.*math.pi # unit rad^2

# load the dataset
nc_base = Dataset(args.input, 'r')
base = projections.DataGrid()
base.y = nc_base.variables['y1']
base.x = nc_base.variables['x1']
base.dy = base.y[1]-base.y[0]
base.dx = base.x[1]-base.x[0]
base.ny = base.y[:].shape[0]
base.nx = base.x[:].shape[0]
base.N = base.ny*base.nx
base.make_grid()

base.lon_grid = nc_base.variables['lon']
base.lat_grid = nc_base.variables['lat']

base.ll_y = base.y_grid.flatten() - base.dy/2.  
base.ll_x = base.x_grid.flatten() - base.dx/2. 

base.lr_y = base.y_grid.flatten() - base.dy/2.  
base.lr_x = base.x_grid.flatten() + base.dx/2. 

base.ur_y = base.y_grid.flatten() + base.dy/2.  
base.ur_x = base.x_grid.flatten() + base.dx/2. 

base.ul_y = base.y_grid.flatten() + base.dy/2.  
base.ul_x = base.x_grid.flatten() - base.dx/2. 

#FIXME: Pick projection based on input file!
base.ll_lon, base.ll_lat = proj_eigen_gl04c(base.ll_x, base.ll_y, inverse=True)
base.lr_lon, base.lr_lat = proj_eigen_gl04c(base.lr_x, base.lr_y, inverse=True)
base.ur_lon, base.ur_lat = proj_eigen_gl04c(base.ur_x, base.ur_y, inverse=True)
base.ul_lon, base.ul_lat = proj_eigen_gl04c(base.ul_x, base.ul_y, inverse=True)

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

nc_scrip = Dataset(lc_scrip, 'w') 
nc_scrip.createDimension('grid_size', base.N)
nc_scrip.createDimension('grid_corners', 4)
nc_scrip.createDimension('grid_rank', 1)
#FIXME: get from input file!
nc_scrip.title = 'CISM EPSG:3413 Grid'
nc_scrip.note = 'Created by Joseph H. Kennedy using pyproj.'

scrip = projections.DataGrid()
scrip.dims = nc_scrip.createVariable('grid_dims','i4', ('grid_rank'))
scrip.dims[:] = base.N

scrip.imask = nc_scrip.createVariable('grid_imask','i4', ('grid_size'))
scrip.imask[:] = numpy.ones(base.N) 
scrip.imask.units = 'unitless'

scrip.center_lat = nc_scrip.createVariable('grid_center_lat','f4', ('grid_size'))
scrip.center_lat[:] = base.lat_grid[0,:,:].ravel()
scrip.center_lat.setncattr('units', base.lat_grid.getncattr('units'))

scrip.center_lon = nc_scrip.createVariable('grid_center_lon','f4', ('grid_size'))
scrip.center_lon[:] = base.lon_grid[0,:,:].ravel()
scrip.center_lon.setncattr('units', base.lon_grid.getncattr('units'))

scrip.corner_lat = nc_scrip.createVariable('grid_corner_lat','f4', ('grid_size','grid_corners',))
scrip.corner_lat[:,:] = base.corner_lat[:,:]
scrip.corner_lat.units = 'degrees'

scrip.corner_lon = nc_scrip.createVariable('grid_corner_lon','f4', ('grid_size','grid_corners',))
scrip.corner_lon[:,:] = base.corner_lon[:,:]
scrip.corner_lon.units = 'degrees'

scrip.area = nc_scrip.createVariable('grid_area','f4', ('grid_size',))
scrip.area[:] = (base.area[:]/surface_area_earth)*sqr_deg_on_sphere
scrip.area.units = 'square degrees'

nc_base.close()
nc_scrip.close()

