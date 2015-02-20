import numpy as np
import pyproj
import datetime
from netCDF4 import Dataset
import scipy
from scipy import interpolate

from ncfunc import copy_atts 

"""
Build a CISM dataset
"""
# get the current time
stamp = datetime.date.today().strftime("%Y_%m_%d")
basename='greenland_500km_'+stamp+'.mcb.nc'

# make some projections
#
#
#

# create base data file
nc_base = Dataset('greenland_500m.mcb.nc','w')


# load in datasets
nc_bamber   = Dataset( 'data/BamberDEM/Greenland_bedrock_topography_V3.nc', 'r') 
nc_seaRise  = Dataset( 'data/SeaRise/Greenland1km.nc', 'r')
nc_racmo2p0 = Dataset( 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc', 'r')


#===== Bamber DEM =====
# this is a 1km dataset
#======================
bamber_y = nc_bamber.variables['projection_y_coordinate']
bamber_ny = bamber_y[:].shape[0] # number of y points for 1km grid

bamber_x = nc_bamber.variables['projection_x_coordinate']
bamber_nx = bamber_x[:].shape[0] # number of x points for 1km grid

# make bamber 500m grid for base
base_y_array = np.arange( min(bamber_y), max(bamber_y)+500., 500., 'f4' )
base_ny = base_y_array.shape[0]

base_x_array = np.arange( min(bamber_x), max(bamber_x)+500., 500., 'f4' )
base_nx = base_x_array.shape[0]

# create our base dimensions
nc_base.createDimension('y',base_ny)
nc_base.createDimension('x',base_nx)

# create some base variables
base_y = nc_base.createVariable('y', 'f4', 'y')
base_y[:] = base_y_array[:]
copy_atts(bamber_y, base_y)

base_x = nc_base.createVariable('x', 'f4', 'x')
base_x[:] = base_x_array[:]
copy_atts(bamber_x, base_x)

# create some grids for interpolation
bamber_y_grid, bamber_x_grid = scipy.meshgrid( bamber_x[:], bamber_y[:] ) #yx then xy correct
base_y_grid, base_x_grid = scipy.meshgrid(base_x[:], base_y[:])           #yx then xy correct

#==== SeaRise Data ====
# this is a 1km dataset
#======================
seaRise_y = nc_seaRise.variables['y']
seaRise_ny = seaRise_y[:].shape[0]

seaRise_x = nc_seaRise.variables['x']
seaRise_nx = seaRise_x[:].shape[0]

# to convert to Bamber 1km grid
seaRise_data = np.ndarray( (bamber_ny, bamber_nx) )
seaRise_y_equal_bamber = 100
seaRise_x_equal_bamber = 500

# to interpolate to Bamber 500m grid
seaRise_500m_data = np.ndarray( (base_ny, base_nx) )

# get basal heat flux
#--------------------
seaRise_data[:,:] = 0.
seaRise_bheatflx = nc_seaRise.variables['bheatflx']
seaRise_data[ seaRise_y_equal_bamber:seaRise_y_equal_bamber+seaRise_ny , seaRise_x_equal_bamber:seaRise_x_equal_bamber+seaRise_nx ] = -seaRise_bheatflx[0,:,:] # invert sign!

# interpolate to 500m
print(seaRise_data.shape)
print(len(bamber_y))
print(len(bamber_x))

seaRise_to_base = interpolate.RectBivariateSpline( bamber_y, bamber_x, seaRise_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
for ii in range(0, base_nx):
    seaRise_500m_data[:,ii] = seaRise_to_base.ev(base_y_grid[:,ii], base_x_grid[:,ii] )

base_bheatflx = nc_base.createVariable('bheatflx', 'f4', ('y','x',) )
base_bheatflx[:,:] = seaRise_500m_data[:,:]
copy_atts(seaRise_bheatflx, base_bheatflx)

# get annual mean air temperature (2m)
#-------------------------------------
seaRise_data[:,:] = 0.
seaRise_presartm = nc_seaRise.variables['presartm']
seaRise_data[ seaRise_y_equal_bamber:seaRise_y_equal_bamber+seaRise_ny , seaRise_x_equal_bamber:seaRise_x_equal_bamber+seaRise_nx ] = seaRise_presartm[0,:,:]

# interpolate to 500m
seaRise_to_base = interpolate.RectBivariateSpline( bamber_y, bamber_x, seaRise_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
for ii in range(0, base_nx):
    seaRise_500m_data[:,ii] = seaRise_to_base.ev(base_y_grid[:,ii], base_x_grid[:,ii] )

base_artm = nc_base.createVariable('artm', 'f4', ('y','x',) )
base_artm[:,:] = seaRise_500m_data[:,:]
copy_atts(seaRise_presartm, base_artm)



#==== RACMO2.0 Data ====
# this is a ?km dataset 
#=======================

