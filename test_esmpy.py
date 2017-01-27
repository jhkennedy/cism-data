#!/usr/bin/env python2

import os
import ESMF
import numpy as np

import argparse
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt

from util import speak
from util import projections
from util.ncfunc import get_nc_file

ESMF.Manager(debug=True)

epsg_cism  = 'greenland_1km_test.nc'
epsg_scrip = 'SCRIPgrid_greenland_1km_test.nc'

lc_base     = 'templates/greenland_1km.epsg3413.nc'
lc_epsg     = 'data/EPSG3413/EPSG3413grid.json'        #NOTE: created by plot_projections.py
lc_epsg_shr = 'data/EPSG3413/EPSG3413grid_shrunk.json' #NOTE: created by plot_grids.py
lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_InSAR    = 'data/InSAR/Joughin2015/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
epsg_InSAR  = 'data/InSAR/Joughin2015/SCRIPgrid_greenland_vel_mosaic500.nc'

def abs_existing_file(file):
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        print("Error! File does not exist: \n    "+file)
        sys.exit(1)
    return file

# parse the command line arguments
parser = argparse.ArgumentParser()   # -h or --help automatically included!

parser.add_argument('-e', '--extended', help='Produce the extended grid.', action='store_true')

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

from data import epsg3413
f_epsg = abs_existing_file(lc_epsg)
speak.verbose(args,"   Found EPSG:3413 grid specs")
f_epsg_shr = abs_existing_file(lc_epsg_shr)
speak.verbose(args,"   Found shrunken EPSG:3413 grid specs")

from data import bamberdem
nc_bamber = get_nc_file(lc_bamber,'r')
speak.verbose(args,"   Found Bamber DEM")
    
from data import insar
try:
    nc_insar = get_nc_file(lc_InSAR,'r')
except Exception:
    speak.verbose(args,"\n   Building InSAR velocity dataset...\n")
    subprocess.call("python util/convert_velocities.py "+os.path.dirname(lc_InSAR), shell=True)
    nc_insar = get_nc_file(lc_InSAR,'r')
speak.verbose(args,"   Found InSAR data")



insar_grid = ESMF.Grid(filename=os.path.join(os.getcwd(),epsg_InSAR), filetype=ESMF.FileFormat.SCRIP)
insar_field = ESMF.Field(insar_grid, staggerloc=ESMF.StaggerLoc.CENTER)
insar_field.data[:,:] = nc_insar.variables['vx'][:,:] 

cism_grid = ESMF.Grid(filename=os.path.join(os.getcwd(),epsg_scrip), filetype=ESMF.FileFormat.SCRIP)
cism_field = ESMF.Field(cism_grid, staggerloc=ESMF.StaggerLoc.CENTER)

print('regridding...')
regrid = ESMF.Regrid(insar_field, cism_field, regrid_method=ESMF.RegridMethod.BILINEAR, 
                     unmapped_action=ESMF.UnmappedAction.IGNORE)
cism_field = regrid(insar_field, cism_field)


print('plotting...')
plt.figure(1, figsize=(10,8), dpi=150)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

nc_base = get_nc_file(lc_base, 'r')
base = projections.DataGrid()
base.y = nc_base.variables['y']
base.x = nc_base.variables['x']
base.ny = base.y[:].shape[0]
base.nx = base.x[:].shape[0] 
base.make_grid()

proj_epsg3413, proj_eigen_gl04c = projections.greenland()

B_ll = ( np.amin(base.x_grid[1,:]), np.amin(base.y_grid[1,:]) )
B_lr = ( np.amax(base.x_grid[1,:]), np.amax(base.y_grid[1,:]) )
B_ur = ( np.amax(base.x_grid[-1,:]), np.amax(base.y_grid[-1,:]) )
B_ul = ( np.amin(base.x_grid[-1,:]), np.amin(base.y_grid[-1,:]) )

B_xs = np.array([ B_ll[0], B_lr[0], B_ur[0], B_ul[0], B_ll[0] ])
B_ys = np.array([ B_ll[1], B_lr[1], B_ur[1], B_ul[1], B_ll[1] ])

B_lons, B_lats = proj_epsg3413(B_xs, B_ys, inverse=True) 
map_ll_lon, map_ll_lat = proj_epsg3413(B_xs[0], B_ys[0], inverse=True) 
map_ur_lon, map_ur_lat = proj_epsg3413(B_xs[2], B_ys[2], inverse=True) 

lon_grid, lat_grid = proj_epsg3413(base.x_grid.ravel(), base.y_grid.ravel(), inverse=True)
lon_grid.shape = base.x_grid.shape 
lat_grid.shape = base.x_grid.shape

glmap = Basemap(llcrnrlon=map_ll_lon, llcrnrlat=map_ll_lat,
                urcrnrlon=map_ur_lon, urcrnrlat=map_ur_lat,
                resolution='l', epsg=3413)

map_x, map_y = glmap(lon_grid, lat_grid)
glmap.pcolormesh(map_x, map_y, cism_field.data)

glmap.fillcontinents(color='gray', lake_color='white')
glmap.drawcoastlines()
glmap.colorbar()


plt.tight_layout()
plt.show()

print('Done!')
