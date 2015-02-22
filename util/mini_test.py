import os
import math
import scipy
import pyproj
import datetime
import subprocess
import numpy as np
from netCDF4 import Dataset
from scipy import interpolate

from ncfunc import copy_atts 

# get the current time
stamp = datetime.date.today().strftime("%Y_%m_%d")
#basename='greenland_500km_'+stamp+'.mcb.nc'

# create base data file
#nc_base = Dataset('play_base.nc','w')
#nc_base = Dataset('play_insar.nc','w')

# load in datasets
nc_bamber   = Dataset( 'data/BamberDEM/Greenland_bedrock_topography_V3.nc', 'r') 
nc_coords   = Dataset( 'util/old/coords.nc', 'r') 

if not ( os.path.exists('data/InSAR/Joughin2012/greenland_vel_mosaic500.nc') ):
    subprocess.call("python util/convert_velocities.py", shell=True)
nc_insar    = Dataset( 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' , 'r')

#==== Projections ====
# All the projections 
# needed for the data 
#=====================
proj_latlong  = pyproj.Proj(proj='latlong', datum='WGS84') 

proj_bamber   = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=-39.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84 +units=m') 
proj_wat      = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=-39.0 +k_0=1.0 +ellps=WGS84 +units=m') 
    # basic bamber projection
proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m')
    # InSAR data in this projections

# EIGEN-GL04C referenced data
#----------------------------
# unfortunately, bed, surface, and thickness data is referenced to EIGEN-GL04C which doesn't exist in proj4. However, EGM2008 should be within ~1m everywhere
# (and within 10-20 cm in most places) so we use the egm08 projection which is available in proj4
if not ( os.path.exists('data/BamberDEM/egm08_25.gtx') ):
    raise Exception("No data/BamberDEM/egm08_25.gtx ! Get it here: http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx") 
proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./data/BamberDEM/egm08_25.gtx')
 
#===== Bamber DEM =====
# this is a 1km dataset
#======================
bamber_y = nc_bamber.variables['projection_y_coordinate']
bamber_ny = bamber_y[:].shape[0] # number of y points for 1km grid

bamber_x = nc_bamber.variables['projection_x_coordinate']
bamber_nx = bamber_x[:].shape[0] # number of x points for 1km grid

# make bamber 500m grid for base
base_ny = bamber_ny
base_nx = bamber_nx

# create some base variables
base_y = np.array( base_nx )
base_y = bamber_y[:]

base_x = np.array( base_nx )
base_x = bamber_x[:]

# create some grids for interpolation
base_y_grid, base_x_grid = scipy.meshgrid(base_y[:], base_x[:], indexing='ij') 


#==== InSAR velocity Data ====
# this is a 500m dataset in   
# the ESPG-3413 projection    
#=============================
insar_y = nc_insar.variables['y']
insar_ny = insar_y[:].shape[0]

insar_x = nc_insar.variables['x']
insar_nx = insar_x[:].shape[0]

insar_vy = nc_insar.variables['vy']

# transform meshes
insar_y_grid, insar_x_grid = scipy.meshgrid( insar_y[:], insar_x[:] , indexing='ij')

insar_data = np.ndarray( (insar_ny, insar_nx) )
insar_data[:,:] = insar_vy[:,:]

#==== Coordinates Transformed ====
# c progams to transform them
#=================================

coords_y_grid = nc_coords.variables['y']
coords_x_grid = nc_coords.variables['x']
coords_ny = insar_ny
coords_nx = insar_nx

#==== recreate trans. ====
# why this no worky???
#=========================

trans_x_grid, trans_y_grid = pyproj.transform(proj_epsg3413, proj_wat, insar_x_grid.flatten(), insar_y_grid.flatten())
trans_y_grid = trans_y_grid.reshape((insar_ny,insar_nx))
trans_x_grid = trans_x_grid.reshape((insar_ny,insar_nx))

trans_y = trans_y_grid[:,0]
trans_x = trans_x_grid[0,:]
coords_y = coords_y_grid[:,0]
coords_x = coords_x_grid[0,:]

#import matplotlib.pyplot as plt 
#plt.imshow( insar_data[::-1,:])
#plt.show()

#base_corners_y = [base_y[0],base_y[0],base_y[-1],base_y[-1] ]
#base_corners_x = [base_x[0],base_x[-1],base_x[0],base_x[-1] ]
#coords_corners_y = [coords_y[0],coords_y[0],coords_y[-1],coords_y[-1] ]
#coords_corners_x = [coords_x[0],coords_x[-1],coords_x[0],coords_x[-1] ]
#trans_corners_y = [trans_y[0],trans_y[0],trans_y[-1],trans_y[-1] ]
#trans_corners_x = [trans_x[0],trans_x[-1],trans_x[0],trans_x[-1] ]
#
#import matplotlib.pyplot as plt
#plt.figure(1)
#plt.scatter( base_corners_y, base_corners_x, c='r' )
#plt.scatter( coords_corners_y, coords_corners_x, s=80, c='b' )
#plt.scatter( trans_corners_y, trans_corners_x, c='g' )
#plt.show()

#======ITWORKS====
# WEHOOOOOOOODLEE 
#=================
# now what the bleep is wrong with my interpolation???



#print(base_corners_y)
#print(base_corners_x)
#
#print(trans_corners_y)
#print(trans_corners_x)

insar_to_base = interpolate.RectBivariateSpline( trans_y[:], trans_x[:], insar_data, s=0) # regular 2d linear interp. but faster

base_data = np.ndarray( (base_ny,base_nx) )
base_data[:,:] = 0.
for ii in range(0, base_nx):
    base_data[:,ii] = insar_to_base.ev(base_y_grid[:,ii], base_x_grid[:,ii] )

import matplotlib.pyplot as plt
plt.figure(2)
plt.imshow(base_data[::-1,:])
plt.show()

base_data = np.ndarray( (base_ny,base_nx) )
base_data[:,:] = 0.
insar_to_base = interpolate.interp2d( trans_y[:], trans_x[:], insar_data) # regular 2d linear interp. but faster
base_data[:,:] = insar_to_base(base_y, base_x )

import matplotlib.pyplot as plt
plt.figure(3)
plt.imshow(base_data[::-1,:])
plt.show()








#nc_base.close()
