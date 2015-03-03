import os
import sys
import datetime
import subprocess
import argparse

import math
import scipy
import pyproj
import numpy as np
from netCDF4 import Dataset
from scipy import interpolate

from util.ncfunc import copy_atts
from util import speak


"""
Build a CISM dataset
"""
#==== Data Locations ====
# Link data here or edit 
#========================
lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_seaRise  = 'data/SeaRise/Greenland1km.nc'
lc_racmo2p0 = 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc'
lc_InSAR    = 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'
lc_mask     = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'


#==== SETUP ====
# get args, time
# load data sets
#===============
stamp = datetime.date.today().strftime("%Y_%m_%d")

# a dummy class for the grid structure
class grid():
    pass

# parse the command line arguments
# -h or --help automatically included!
parser = argparse.ArgumentParser()

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

speak.notquiet(args,"\nBuilding the Greenland datasets in the Bamber projection.")
speak.notquiet(args,  "=========================================================\n")

# create base data file
f_base = 'templates/greenland_1km.mcb.nc'
nc_base = Dataset(f_base,'w')
speak.verbose(args,"Building the base dataset: "+f_base+"\n")

# load in datasets
speak.notquiet(args,"Loading the datasets.")

nc_bamber   = Dataset( lc_bamber, 'r') 
from data import bamberdem
speak.verbose(args,"   Found Bamber DEM")

nc_seaRise  = Dataset( lc_seaRise, 'r')
from data import searise
speak.verbose(args,"   Found SeaRise data")

nc_racmo2p0 = Dataset( lc_racmo2p0, 'r')
from data import racmo2p0
speak.verbose(args,"   Found RACMO 2.0 data")

if not ( os.path.exists(lc_InSAR) ):
    speak.verbose(args,"\n   Building InSAR velocity dataset...\n")
    subprocess.call("python util/convert_velocities.py "+os.path.dirname(lc_InSAR), shell=True)
nc_insar    = Dataset( lc_InSAR , 'r')
from data import insar
speak.verbose(args,"   Found InSAR data")

nc_massCon  = Dataset(lc_massCon,'r')
from data import icebridge
speak.verbose(args,"   Found Mass Conserving Bed data")

nc_mask     = Dataset( lc_mask, 'r'  )
from data import ice2sea
speak.verbose(args,"   Found Zurich mask")

speak.verbose(args,"\n   All data files found!")
#===== Bamber DEM =====
# this is a 1km dataset
#======================
speak.notquiet(args,"\nBuilding the base grid."),

base = grid()
base = bamberdem.build_base(nc_bamber, nc_base, base)

speak.notquiet(args,"   Done!")

#==== Projections ====
# All the projections 
# needed for the data 
#=====================
speak.notquiet(args,"\nBuilding the projections.")
    
proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m')
    # InSAR data in this projections

# EIGEN-GL04C referenced data:
#----------------------------
# unfortunately, bed, surface, and thickness data is referenced to EIGEN-GL04C which doesn't exist in proj4. However, EGM2008 should be within ~1m everywhere
# (and within 10-20 cm in most places) so we use the egm08 projection which is available in proj4
path_bamber = os.path.dirname(lc_bamber)
if not ( os.path.exists(path_bamber+'/egm08_25.gtx') ):
    raise Exception("No "+path_bamber+"/egm08_25.gtx ! Get it here: http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx") 
#NOTE: Bamber projection appears to not actually have any fasle northings or eastings. 
#proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')
proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')

# transform meshes. 
speak.verbose(args,"   Creating the transform meshes: base Bamber grid to EPSG-3413.")

trans = grid()
trans.ny = base.ny
trans.nx = base.nx

trans.x_grid, trans.y_grid = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, base.x_grid.flatten(), base.y_grid.flatten())
trans.y_grid = trans.y_grid.reshape((base.ny,base.nx))
trans.x_grid = trans.x_grid.reshape((base.ny,base.nx))

speak.notquiet(args,"   Done!")
#==== SeaRise Data ====
# this is a 1km dataset
#======================
speak.notquiet(args,"\nGetting bheatflx and presartm from the SeaRise data.")
searise.get_bheatflx_artm(args, nc_seaRise, nc_base, base)

#==== RACMO2.0 Data ====
# this is a 1km dataset 
#=======================
speak.notquiet(args,"\nGetting acab from the RACMO 2.0 data.")
racmo2p0.get_acab(args, nc_racmo2p0, nc_base, base)

#==== InSAR velocity Data ====
# this is a 500m dataset in   
# the ESPG-3413 projection    
#=============================
speak.notquiet(args,"\nGetting vy, vx, ey, and ex from the InSAR data.")
insar.get_velocity(args, nc_insar, nc_base, trans)

#==== Mass Conserving Bed Data ===
# This is the new (2015) bed data 
#=================================
# new vars: 'thk' 'topg' 'topgerr' 'usrf'** make from thk+topg?? 'usrfRMSE'** doesn't actually exist.
speak.notquiet(args,"\nGetting thk, topg, and topgerr from the mass conserving bed data.")
icebridge.get_mcb(args, nc_massCon, nc_base, trans, proj_eigen_gl04c, proj_epsg3413)

nc_base.close()
#==== Zurich mask ====
# apply mask, and get 
# new surface variable
#=====================
speak.notquiet(args,"\nGetting the Zurich Mask.")

nc_base = Dataset(f_base,'r+')
base_thk  = nc_base.variables['thk']
thk_data = base_thk[:,:]

base_topg = nc_base.variables['topg']
topg_data = base_topg[:,:]

mask = nc_mask.variables['IceSheetMask']

speak.verbose(args,"   Applying mask to thk.")
thk_data = thk_data * mask[:,:]
base_thk[:,:] = thk_data

speak.verbose(args,"   Creating usrf.")
base_usrf = nc_base.createVariable('usrf', 'f4',('y','x',))
base_usrf[:,:] = thk_data + topg_data

nc_mask.close()
nc_base.close()

#==== add time dim and shrink ====
# apply to all the variables and  
# shrink to size around ice sheet 
#=================================
speak.notquiet(args,"\nAdding the time dimension and creating the 1km dataset.")

# shrink dataset to the ice sheet
y_shrink = [100,2900+1] #NOTE: python stop exclusive, nco stop inclusive!
x_shrink = [500,2000+1] #NOTE: python stop exclusive, nco stop inclusive!

nc_base = Dataset(f_base,'r')
base.y = nc_base.variables['y']
base.ny = base.y[ y_shrink[0]:y_shrink[1] ].shape[0]
base.x = nc_base.variables['x']
base.nx = base.x[ x_shrink[0]:x_shrink[1] ].shape[0]

f_1km = 'complete/greenland_1km_'+stamp+'.mcb.nc'
speak.verbose(args,"   Writing "+f_1km)
nc_1km = Dataset(f_1km, 'w')
nc_1km.createDimension('time', None )
nc_1km.createDimension('y1', base.ny)
nc_1km.createDimension('x1', base.nx)

time = nc_1km.createVariable('time', 'f4', ('time',))
y1   = nc_1km.createVariable('y1',   'f4', ('y1',)  )
x1   = nc_1km.createVariable('x1',   'f4', ('x1',)  )

copy_atts(base.y, y1)
copy_atts(base.x, x1)

y1[:] = base.y[ y_shrink[0]:y_shrink[1] ]
x1[:] = base.x[ x_shrink[0]:x_shrink[1] ]
time[0] = 0.

for var_name, var_data in nc_base.variables.iteritems() : 
    if (var_name != 'x' and var_name != 'y'):
        var_1km = nc_1km.createVariable(var_name, 'f4', ('time','y1','x1',))
        var_1km[0,:,:] = var_data[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]
        copy_atts(var_data, var_1km)
#copy_atts(var_data, var_1km)

speak.verbose(args,"   Writing the 1km config file.")
# write 1km config file
config_dict= {'REPLACE_EWN':str(base.nx), 
              'REPLACE_NSN':str(base.ny), 
              'REPLACE_DEW':str(1000), 
              'REPLACE_DNS':str(1000), 
              'REPLACE_NAME':'greenland_1km_'+stamp+'.mcb.nc', 
              'REPLACE_OUT':'greenland_1km_'+stamp+'.mcb.out.nc', 
              'REPLACE_KM':'1 km' }

base_config = open('templates/greenland_base.mcb.config','r')
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
# make 2, 4 and 8 km dataset
#===========================
speak.notquiet(args,"\nCreating coarser datasets.")

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
    
    f_coarse = 'complete/greenland_'+coarse_list[ii]+'_'+stamp+'.mcb.nc'
    speak.verbose(args,"   Writing "+f_coarse)
    nc_coarse = Dataset( f_coarse,'w' )
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
    # set file permissions
    os.chmod(f_coarse, 0o644) # uses an Octal number!

    # write config files
    speak.verbose(args,"   Writing the "+coarse_list[ii]+" config file.")
    config_dict= {'REPLACE_EWN':str(nx_coarse), 
                  'REPLACE_NSN':str(ny_coarse), 
                  'REPLACE_DEW':str(skip*1000), 
                  'REPLACE_DNS':str(skip*1000), 
                  'REPLACE_NAME':'greenland_'+coarse_list[ii]+'_'+stamp+'.mcb.nc', 
                  'REPLACE_OUT':'greenland_'+coarse_list[ii]+'_'+stamp+'.mcb.out.nc', 
                  'REPLACE_KM':str(skip)+' km' }
    
    base_config = open('templates/greenland_base.mcb.config','r')
    
    f_config = 'complete/greenland_'+coarse_list[ii]+'.mcb.config'
    out_config  = open(f_config,'w')
    for line in base_config :
        for src, target in config_dict.iteritems() :
            line = line.replace(src, target)
        out_config.write(line)
    
    base_config.close()
    out_config.close()
    # set file permissions
    os.chmod(f_config, 0o755) # uses an octal number!


nc_1km.close()

#==== and done! ====
speak.notquiet(args,"\nFinished building the datasets.")
