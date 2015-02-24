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

# create base data file
nc_base = Dataset('greenland_1km.mcb.nc','w')

# load in datasets
nc_bamber   = Dataset( 'data/BamberDEM/Greenland_bedrock_topography_V3.nc', 'r') 
nc_seaRise  = Dataset( 'data/SeaRise/Greenland1km.nc', 'r')
nc_racmo2p0 = Dataset( 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc', 'r')

if not ( os.path.exists('data/InSAR/Joughin2012/greenland_vel_mosaic500.nc') ):
    subprocess.call("python util/convert_velocities.py", shell=True)
nc_insar    = Dataset( 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' , 'r')

nc_massCon  = Dataset('data/IceBridge/Greenland/MCdataset-2014-11-19.nc','r')
nc_mask     = Dataset( 'data/Ice2Sea/greenland_geometry_icesheet_mask_Zurich.nc' , 'r'  )


#==== Projections ====
# All the projections 
# needed for the data 
#=====================
#NOTE: Bamber projection appears to not actually have any fasle northings or eastings. 
proj_latlong  = pyproj.Proj(proj='latlong', datum='WGS84') 
    # needed to convert between projections
#proj_bamber   = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=-39.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +ellps=WGS84 +units=m') 
proj_bamber   = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=-39.0 +k_0=1.0 +ellps=WGS84 +units=m') 
    # basic bamber projection
proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m')
    # InSAR data in this projections

# EIGEN-GL04C referenced data
#----------------------------
# unfortunately, bed, surface, and thickness data is referenced to EIGEN-GL04C which doesn't exist in proj4. However, EGM2008 should be within ~1m everywhere
# (and within 10-20 cm in most places) so we use the egm08 projection which is available in proj4
if not ( os.path.exists('data/BamberDEM/egm08_25.gtx') ):
    raise Exception("No data/BamberDEM/egm08_25.gtx ! Get it here: http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx") 
#proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./data/BamberDEM/egm08_25.gtx')
proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +geoidgrids=./data/BamberDEM/egm08_25.gtx')

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
base_y_grid, base_x_grid = scipy.meshgrid(base_y[:], base_x[:], indexing='ij') 

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

# transform meshes
insar_y_grid, insar_x_grid = scipy.meshgrid( insar_y[:], insar_x[:] , indexing='ij')

trans_x_grid, trans_y_grid = pyproj.transform(proj_epsg3413, proj_eigen_gl04c, insar_x_grid.flatten(), insar_y_grid.flatten())
trans_y_grid = trans_y_grid.reshape((insar_ny,insar_nx))
trans_x_grid = trans_x_grid.reshape((insar_ny,insar_nx))

trans_y = trans_y_grid[:,0]
trans_x = trans_x_grid[0,:]

insar_data = np.ndarray( (insar_ny, insar_nx) )
base_data = np.ndarray( (base_ny,base_nx) )


for vv in ['vy','vx','ey','ex'] :
    insar_data[:,:] = 0.
    base_data[:,:] = 0.
    
    insar_var = nc_insar.variables[ vv ]
    insar_data[:,:] = insar_var[:,:]

    insar_to_base = interpolate.RectBivariateSpline( trans_y[:], trans_x[:], insar_data, s=0) # regular 2d linear interp. but faster

    for ii in range(0, base_nx):
        base_data[:,ii] = insar_to_base.ev(base_y_grid[:,ii], base_x_grid[:,ii] )

    base_var = nc_base.createVariable( vv, 'f4', ('y','x',) )
    base_var[:,:] = base_data[:,:]
    copy_atts(insar_var, base_var)

# drop temp. variables
trans_y_grid = None
trans_x_grid = None
trans_y = None
trans_x = None
nc_insar.close()
#==== Mass Conserving Bed Data ===
# This is the new (2015) bed data 
#=================================
# new vars: 'thk' 'topg' 'topgerr' 'usrf'** make from thk+topg?? 'usrfRMSE'** doesn't actually exist.
massCon_y = nc_massCon.variables['y']
massCon_ny = massCon_y[:].shape[0]

massCon_x = nc_massCon.variables['x']
massCon_nx = massCon_x[:].shape[0]

# transform meshes
massCon_y_grid, massCon_x_grid = scipy.meshgrid(massCon_y[::-1], massCon_x[:], indexing='ij')  # y fliped when compaired to Bamber

#massCon_bed       = nc_massCon.variables['bed']
#massCon_thickness = nc_massCon.variables['thickness']
#massCon_surface   = nc_massCon.variables['surface']
#massCon_errbed    = nc_massCon.variables['errbed']
#massCon_mask      = mc_massCon.variables['mask']

massCon_data = np.ndarray( (massCon_ny,massCon_nx) )
trans_data   = np.ndarray( (massCon_ny,massCon_nx) )
base_data    = np.ndarray( (base_ny,base_nx) )

var_list       = [ 'thickness', 'bed' ]
rename_massCon = [ 'thk',       'topg']
for vv in range(0, len(var_list) ) :
    var = var_list[vv]
    massCon_data[:,:] = 0.
    trans_data[:,:]   = 0.
    base_data[:,:]    = 0.

    massCon_var = nc_massCon.variables[var]
    massCon_data[:,:] = massCon_var[::-1,:] # y fliped when compaired to Bamber

    trans_x_grid, trans_y_grid, trans_data = pyproj.transform(proj_epsg3413, proj_eigen_gl04c, massCon_x_grid.flatten(), massCon_y_grid.flatten(), massCon_data.flatten())
    trans_data   =   trans_data.reshape( (massCon_ny,massCon_nx) )
    trans_y_grid = trans_y_grid.reshape( (massCon_ny,massCon_nx) ) 
    trans_x_grid = trans_x_grid.reshape( (massCon_ny,massCon_nx) ) 

    trans_y = trans_y_grid[:,0]
    trans_x = trans_x_grid[0,:]
    

    massCon_to_base = interpolate.RectBivariateSpline( trans_y[:], trans_x[:], trans_data, s=0) # regular 2d linear interp. but faster

    for ii in range(0, base_nx):
        base_data[:,ii] = massCon_to_base.ev(base_y_grid[:,ii], base_x_grid[:,ii] )
    
    
    base_var = nc_base.createVariable( rename_massCon[vv], 'f4', ('y','x',) )
    base_var[:,:] = base_data[:,:]
    #copy_atts(massCon_var, base_var)
    #NOTE: this data is formated short, while the rest of the data is a float. Does not like _FillValue from short data. 
    atts = massCon_var.ncattrs()
    for ii in range(len(atts)):
        if (atts[ii] != '_FillValue'):
            base_var.setncattr(atts[ii], massCon_var.getncattr(atts[ii]))
        #else:
        #    base_var.setncattr('_FillValue', -9999.)

# Now get bed error -- not in for loop because not z data!
#---------------------------------------------------------
base_data[:,:]    = 0.

massCon_errbed = nc_massCon.variables['errbed']

#trans_x_grid, trans_y_grid = pyproj.transform(proj_epsg3413, proj_eigen_gl04c, massCon_x_grid.flatten(), massCon_y_grid.flatten())
#trans_y_grid = trans_y_grid.reshape( (massCon_ny,massCon_nx) ) 
#trans_x_grid = trans_x_grid.reshape( (massCon_ny,massCon_nx) ) 
#
#trans_y = trans_y_grid[:,0]
#trans_x = trans_x_grid[0,:]

massCon_to_base = interpolate.RectBivariateSpline( trans_y[:], trans_x[:], massCon_errbed[::-1,:], s=0) # regular 2d linear interp. but faster
                                                                            # y fliped when compaired to Bamber
for ii in range(0, base_nx):
    base_data[:,ii] = massCon_to_base.ev(base_y_grid[:,ii], base_x_grid[:,ii] )

base_var = nc_base.createVariable( 'topgerr', 'f4', ('y','x',) )
base_var[:,:] = base_data[:,:]
#copy_atts(massCon_errbed, base_var)
#NOTE: this data is formated short, while the rest of the data is a float. Does not like _FillValue from short data. 
atts = massCon_var.ncattrs()
for ii in range(len(atts)):
    if (atts[ii] != '_FillValue'):
        base_var.setncattr(atts[ii], massCon_var.getncattr(atts[ii]))
   # else:
   #     base_var.setncattr('_FillValue', -9999.)

# drop temp. variables
trans_data = None
trans_y_grid = None
trans_x_grid = None
trans_y = None
trans_x = None
nc_massCon.close()


#==== Zurich mask ====
# apply mask, and get 
# new surface variable
#=====================
nc_base.close()

nc_base = Dataset('greenland_1km.mcb.nc','r+')
base_thk  = nc_base.variables['thk']
thk_data = base_thk[:,:]

base_topg = nc_base.variables['topg']
topg_data = base_topg[:,:]

mask = nc_mask.variables['IceSheetMask']

thk_data = thk_data * mask[:,:]
base_thk[:,:] = thk_data

base_usrf = nc_base.createVariable('usrf', 'f4',('y','x',))
base_usrf[:,:] = thk_data + topg_data

nc_mask.close()
nc_base.close()
#==== add time dim ====
# apply to all the vars
#======================
nc_base = Dataset('greenland_1km.mcb.nc','r')
base_y = nc_base.variables['y']
base_ny = base_y[:].shape[0]
base_x = nc_base.variables['x']
base_nx = base_x[:].shape[0]


nc_1km = Dataset('complete/greenland_1km_'+stamp+'.mcb.nc', 'w')
nc_1km.createDimension('time', None )
nc_1km.createDimension('y1', base_ny)
nc_1km.createDimension('x1', base_nx)

time = nc_1km.createVariable('time', 'f4', ('time',))
y1   = nc_1km.createVariable('y1',   'f4', ('y1',)  )
x1   = nc_1km.createVariable('x1',   'f4', ('x1',)  )

copy_atts(base_y, y1)
copy_atts(base_x, x1)

y1[:] = base_y[:]
x1[:] = base_x[:]
time[0] = 0.

for var_name, var_data in nc_base.variables.iteritems() : 
    if (var_name != 'x' and var_name != 'y'):
        var_1km = nc_1km.createVariable(var_name, 'f4', ('time','y1','x1',))
        var_1km[0,:,:] = var_data[:,:]
        copy_atts(var_data, var_1km)
#copy_atts(var_data, var_1km)

# write 1km config file
config_dict= {'REPLACE_EWN':str(base_nx), 
              'REPLACE_NSN':str(base_ny), 
              'REPLACE_DEW':str(1000), 
              'REPLACE_DNS':str(1000), 
              'REPLACE_NAME':'greenland_1km_'+stamp+'.mcb.nc', 
              'REPLACE_OUT':'greenland_1km_'+stamp+'.mcb.out.nc', 
              'REPLACE_KM':'1 km' }

base_config = open('greenland_base.mcb.config','r')
out_config  = open('complete/greenland_1km.mcb.config','w')
for line in base_config :
    for src, target in config_dict.iteritems() :
        line = line.replace(src, target)
    out_config.write(line)

base_config.close()
out_config.close()

nc_1km.close()
nc_base.close()
#==== Stamp and coarsen ====
# make 1, 2, 4 and 8 km data
#===========================
nc_1km = Dataset('complete/greenland_1km_'+stamp+'.mcb.nc','r' )

y_1km = nc_1km.variables['y1']
ny_1km = y_1km[:].shape[0]

x_1km = nc_1km.variables['x1']
nx_1km = x_1km[:].shape[0]


coarse_list = ['2km','4km','8km']
for ii in range(0, len(coarse_list)):
    skip = 2**(ii+1)
    ny_coarse = y_1km[::skip].shape[0]
    nx_coarse = x_1km[::skip].shape[0]

    nc_coarse = Dataset( 'complete/greenland_'+coarse_list[ii]+'_'+stamp+'.mcb.nc','w' )
    nc_coarse.createDimension('time', None )
    nc_coarse.createDimension('y1', ny_coarse)
    nc_coarse.createDimension('x1', nx_coarse)

    time_coarse = nc_coarse.createVariable('time', 'f4', ('time',))
    y1_coarse   = nc_coarse.createVariable('y1',   'f4', ('y1',)  )
    x1_coarse   = nc_coarse.createVariable('x1',   'f4', ('x1',)  )

    copy_atts(y_1km, y1_coarse)
    copy_atts(x_1km, x1_coarse)

    y1_coarse[:] = y_1km[::skip]
    x1_coarse[:] = x_1km[::skip]
    time_coarse[0] = 0.

    for var_name, var_data in nc_1km.variables.iteritems() : 
        if (var_name != 'x1' and var_name != 'y1' and var_name != 'time'):
            var_coarse = nc_coarse.createVariable(var_name, 'f4', ('time','y1','x1',))
            var_coarse[0,:,:] = var_data[0,::skip,::skip]
            copy_atts(var_data, var_coarse)

    nc_coarse.close()

    # write config files
    config_dict= {'REPLACE_EWN':str(nx_coarse), 
                  'REPLACE_NSN':str(ny_coarse), 
                  'REPLACE_DEW':str(skip*1000), 
                  'REPLACE_DNS':str(skip*1000), 
                  'REPLACE_NAME':'greenland_'+coarse_list[ii]+'_'+stamp+'.mcb.nc', 
                  'REPLACE_OUT':'greenland_'+coarse_list[ii]+'_'+stamp+'.mcb.out.nc', 
                  'REPLACE_KM':str(skip)+' km' }

    base_config = open('greenland_base.mcb.config','r')
    out_config  = open('complete/greenland_'+coarse_list[ii]+'.mcb.config','w')
    for line in base_config :
        for src, target in config_dict.iteritems() :
            line = line.replace(src, target)
        out_config.write(line)
    
    base_config.close()
    out_config.close()

nc_1km.close()

