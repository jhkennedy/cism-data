#!/usr/bin/env python

import os
import sys
import math
import numpy
import scipy
import argparse

from netCDF4 import Dataset
from shapely.geometry import shape 

from util import projections

"""
Create grids useful for makeing timeseries.  
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
        help='netCDF dataset to build a climo grid for.')

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

#FIXME: Pick projection based on input file!
lon_grid, lat_grid = proj_eigen_gl04c(base.x_grid.ravel(), base.y_grid.ravel(), inverse=True)
lon_grid.shape = base.x_grid.shape 
lat_grid.shape = base.x_grid.shape

base.lon_grid = lon_grid
base.lat_grid = lat_grid

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

sys.stdout.write("\r   [%-60s] %d%% \n" % ('='*60, 100.))
sys.stdout.flush()
base.area.shape = base.x_grid.shape



base.y0 = (base.y[:-1] + base.y[1:])/2.
base.x0 = (base.x[:-1] + base.x[1:])/2.
base.y0_grid, base.x0_grid = scipy.meshgrid(base.y0[:], base.x0[:], indexing='ij')
base.dim0 = (base.ny-1, base.nx-1)
base.N0 = base.dim0[0]*base.dim0[1]

#FIXME: Pick projection based on input file!
lon0_grid, lat0_grid = proj_eigen_gl04c(base.x0_grid.ravel(), base.y0_grid.ravel(), inverse=True)
lon0_grid.shape = base.x0_grid.shape 
lat0_grid.shape = base.x0_grid.shape

base.lon0_grid = lon0_grid
base.lat0_grid = lat0_grid

base.ll_y0 = base.y0_grid.flatten() - base.dy/2.  
base.ll_x0 = base.x0_grid.flatten() - base.dx/2. 

base.lr_y0 = base.y0_grid.flatten() - base.dy/2.  
base.lr_x0 = base.x0_grid.flatten() + base.dx/2. 

base.ur_y0 = base.y0_grid.flatten() + base.dy/2.  
base.ur_x0 = base.x0_grid.flatten() + base.dx/2. 

base.ul_y0 = base.y0_grid.flatten() + base.dy/2.  
base.ul_x0 = base.x0_grid.flatten() - base.dx/2. 

#FIXME: Pick projection based on input file!
base.ll_lon0, base.ll_lat0 = proj_eigen_gl04c(base.ll_x0, base.ll_y0, inverse=True)
base.lr_lon0, base.lr_lat0 = proj_eigen_gl04c(base.lr_x0, base.lr_y0, inverse=True)
base.ur_lon0, base.ur_lat0 = proj_eigen_gl04c(base.ur_x0, base.ur_y0, inverse=True)
base.ul_lon0, base.ul_lat0 = proj_eigen_gl04c(base.ul_x0, base.ul_y0, inverse=True)

base.corner_lat0 = numpy.column_stack((base.ll_lat0, base.lr_lat0, base.ur_lat0, base.ul_lat0))
base.corner_lon0 = numpy.column_stack((base.ll_lon0, base.lr_lon0, base.ur_lon0, base.ul_lon0))

min_lat0 = numpy.amin(base.corner_lat0)
max_lat0 = numpy.amax(base.corner_lat0)

min_lon0 = numpy.amin(base.corner_lon0)
max_lon0 = numpy.amax(base.corner_lon0)

proj_aea0 = projections.equal_area(min_lat0, max_lat0, (max_lon0+min_lon0)/2.)

# get the area for each grid cell
sys.stdout.write("   [%-60s] %d%%" % ('='*0, 0.))
sys.stdout.flush()
base.area0 = numpy.zeros(base.N0)
for ii in range(base.N0):
    ctr = (ii*60)/base.N0
    if not (ii % 100): 
        sys.stdout.write("\r   [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
        sys.stdout.flush()
    
    lat0 = base.corner_lat0[ii,:]
    lon0 = base.corner_lon0[ii,:]
    x0, y0 = proj_aea(lon0, lat0)
   
    points0 = {'type':'polygon', 'coordinates':[zip(x0,y0)]}
    base.area0[ii] = shape(points0).area #m^2

sys.stdout.write("\r   [%-60s] %d%% \n" % ('='*60, 100.))
sys.stdout.flush()
base.area0.shape = base.x0_grid.shape

path_climo, name_climo = os.path.split(args.input)
lc_climo = os.path.join(path_climo, 'climo_'+name_climo)

nc_climo = Dataset(lc_climo, 'w', format='NETCDF4')
nc_climo.createDimension('time', None )
nc_climo.createDimension('y1', base.ny)
nc_climo.createDimension('x1', base.nx)
nc_climo.createDimension('y0', base.ny-1)
nc_climo.createDimension('x0', base.nx-1)

climo = projections.DataGrid()

climo.time = nc_climo.createVariable('time', 'f4', ('time',))
climo.time.long_name = "Model time"
climo.time.standard_name = "time"
climo.time.units = "common_year since 1-1-1 0:0:0"
climo.time.calendar = "none"
climo.time[0] = 0.

climo.y   = nc_climo.createVariable('y1',   'f4', ('y1',)  )
climo.y.long_name = "Cartesian y-coordinate"
climo.y.units = "meter"
climo.y.axis = "Y"
climo.y[:] = base.y[:]

climo.x   = nc_climo.createVariable('x1',   'f4', ('x1',)  )
climo.x.long_name = "Cartesian x-coordinate, velocity grid"
climo.x.units = "meter"
climo.x.axis = "X"
climo.x[:] = base.x[:]

climo.y0   = nc_climo.createVariable('y0',   'f4', ('y0',)  )
climo.y0.long_name = "Cartesian y-coordinate, velocity grid"
climo.y0.units = "meter"
climo.y0.axis = "Y"
climo.y0[:] = base.y0[:]

climo.x0   = nc_climo.createVariable('x0',   'f4', ('x0',)  )
climo.x0.long_name = "Cartesian x-coordinate, velocity grid"
climo.x0.units = "meter"
climo.x0.axis = "X"
climo.x0[:] = base.x0[:]

climo.lon_grid = nc_climo.createVariable('lon', 'f4', ('time','y1','x1',))
climo.lon_grid.units = 'degrees'
climo.lon_grid.long_name = 'grid center longitude'
climo.lon_grid.note = 'Created by Joseph H. Kennedy using pyproj.'
climo.lon_grid[0,:,:] = base.lon_grid[:,:]

climo.lat_grid = nc_climo.createVariable('lat', 'f4', ('time','y1','x1',))
climo.lat_grid.units = 'degrees'
climo.lat_grid.lat_name = 'grid center latitude'
climo.lat_grid.note = 'Created by Joseph H. Kennedy using pyproj.'
climo.lat_grid[0,:,:] = base.lat_grid[:,:]

climo.lon0_grid = nc_climo.createVariable('lon0', 'f4', ('time','y0','x0',))
climo.lon0_grid.units = 'degrees'
climo.lon0_grid.long_name = 'grid center longitude, velocity grid'
climo.lon0_grid.note = 'Created by Joseph H. Kennedy using pyproj.'
climo.lon0_grid[0,:,:] = base.lon0_grid[:,:]

climo.lat0_grid = nc_climo.createVariable('lat0', 'f4', ('time','y0','x0',))
climo.lat0_grid.units = 'degrees'
climo.lat0_grid.lat_name = 'grid center latitude, velocity grid'
climo.lat0_grid.note = 'Created by Joseph H. Kennedy using pyproj.'
climo.lat0_grid[0,:,:] = base.lat0_grid[:,:]

climo.area = nc_climo.createVariable('area','f4', ('time','y1','x1',))
climo.area[0,:,:] = (base.area[:,:]/surface_area_earth)*sqr_deg_on_sphere
climo.area.units = 'square degrees'

climo.area0 = nc_climo.createVariable('area0','f4', ('time','y0','x0',))
climo.area0[0,:,:] = (base.area0[:,:]/surface_area_earth)*sqr_deg_on_sphere
climo.area0.units = 'square degrees'

nc_base.close()
nc_climo.close()
os.chmod(lc_climo, 0o644)   # uses an octal number!

