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

"""
Build a CISM dataset
"""
# get the current time
stamp = datetime.date.today().strftime("%Y_%m_%d")
#basename='greenland_500km_'+stamp+'.mcb.nc'

# create base data file
nc_base = Dataset('greenland_1km.mcb.nc','w')

# load in datasets
nc_bamber   = Dataset( 'data/BamberDEM/Greenland_bedrock_topography_V3.nc', 'r') 
nc_seaRise  = Dataset( 'data/SeaRise/Greenland1km.nc', 'r')
nc_racmo2p0 = Dataset( 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc', 'r')

if not ( os.path.exists('data/InSAR/Joughin2012/greenland_vel_mosaic500.nc') ):
    subprocess.call("python util/convert_velocities.py", shell=True)
nc_insar    = Dataset( 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' , 'r')

#==== Projections ====
# All the projections 
# needed for the data 
#=====================
proj_latlong  = pyproj.Proj(proj='latlong', datum='WGS84') 
    # needed to convert between projections
proj_bamber   = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84') 
    # basic bamber projection
proj_epsg3413 = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
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

# create our base dimensions
nc_base.createDimension('y',base_ny)
nc_base.createDimension('x',base_nx)

# create some base variables
base_y = nc_base.createVariable('y', 'f4', 'y')
base_y[:] = bamber_y[:]
copy_atts(bamber_y, base_y) #FIXME: units say km, but it's actuall in m

base_x = nc_base.createVariable('x', 'f4', 'x')
base_x[:] = bamber_x[:]
copy_atts(bamber_x, base_x) #FIXME: units say km, but it's actuall in m

# create some grids for interpolation
base_y_grid, base_x_grid = scipy.meshgrid(base_x[:], base_y[:]) #yx then xy is correct

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

# get basal heat flux
#--------------------
seaRise_data[:,:] = 0.
seaRise_bheatflx = nc_seaRise.variables['bheatflx']
seaRise_data[ seaRise_y_equal_bamber:seaRise_y_equal_bamber+seaRise_ny , seaRise_x_equal_bamber:seaRise_x_equal_bamber+seaRise_nx ] = -seaRise_bheatflx[0,:,:] # invert sign!

base_bheatflx = nc_base.createVariable('bheatflx', 'f4', ('y','x',) )
base_bheatflx[:,:] = seaRise_data[:,:] # base = bamber 1km
copy_atts(seaRise_bheatflx, base_bheatflx)

# get annual mean air temperature (2m)
#-------------------------------------
seaRise_data[:,:] = 0.
seaRise_presartm = nc_seaRise.variables['presartm']
seaRise_data[ seaRise_y_equal_bamber:seaRise_y_equal_bamber+seaRise_ny , seaRise_x_equal_bamber:seaRise_x_equal_bamber+seaRise_nx ] = seaRise_presartm[0,:,:]

base_artm = nc_base.createVariable('artm', 'f4', ('y','x',) )
base_artm[:,:] = seaRise_data[:,:] # base = bamber 1km
copy_atts(seaRise_presartm, base_artm)

nc_seaRise.close()
#==== RACMO2.0 Data ====
# this is a 1km dataset 
#=======================
racmo2p0_data = np.ndarray( (base_ny, base_nx) )

racmo2p0_data[:,:] = 0.
racmo2p0_smb = nc_racmo2p0.variables['smb']
racmo2p0_data[:,:] = racmo2p0_smb[:,::-1].transpose() / 910.
racmo2p0_data = np.ma.masked_invalid(racmo2p0_data) # find invalid data and create a mask
racmo2p0_data = racmo2p0_data.filled(0.)            # fill invalid data with zeros

base_acab = nc_base.createVariable( 'acab', 'f4', ('y','x',) )
base_acab[:,:] = racmo2p0_data[:,:]
copy_atts(racmo2p0_smb, base_acab) #FIXME: check atribute units -- divided by 910 earlier

nc_racmo2p0.close()
#==== InSAR velocity Data ====
# this is a 500m dataset in   
# the ESPG-3413 projection    
#=============================
insar_y = nc_insar.variables['y']
insar_ny = insar_y[:].shape[0]

insar_x = nc_insar.variables['x']
insar_nx = insar_x[:].shape[0]

insar_vy = nc_insar.variables['vy']
insar_ey = nc_insar.variables['ey']

insar_vx = nc_insar.variables['vx']
insar_ex = nc_insar.variables['ex']


# transform meshes
insar_y_grid, insar_x_grid = scipy.meshgrid( insar_x[:], insar_y[:] ) #yx then xy correct

insar_data = np.ndarray( (insar_ny, insar_nx) )
insar_data[:,:] = insar_vy[:,:]

trans_y_grid, trans_x_grid = pyproj.transform(proj_bamber, proj_epsg3413, base_y_grid.flatten(), base_x_grid.flatten())
trans_y_grid = trans_y_grid.reshape((base_ny,base_nx))
trans_x_grid = trans_x_grid.reshape((base_ny,base_nx))

####PROBLEMISWITHPROJECTION####


base_to_insar = interpolate.RectBivariateSpline( insar_y[:], insar_x[:], insar_data, s=0) # regular 2d linear interp. but faster

base_data = np.ndarray( (base_ny,base_nx) )
base_data[:,:] = 0.
for ii in range(0, base_nx):
    base_data[:,ii] = base_to_insar.ev(trans_y_grid[:,ii], trans_x_grid[:,ii] )

print(trans_y_grid[1120,1754])
print(trans_x_grid[1120,1754])
print(base_data[1120,1754])
print(base_to_insar(1120,1754))

base_vy = nc_base.createVariable( 'vy', 'f4', ('y','x',) )
base_vy[:,:] = base_data[:,:]
copy_atts(insar_vy, base_vy)
