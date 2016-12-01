#!/usr/bin/env python

import os
import sys
import datetime
import subprocess
import argparse

from util import speak
from util import projections
from util.ncfunc import get_nc_file

"""
Build a CISM dataset
"""

def abs_existing_file(file):
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        print("Error! File does not exist: \n    "+file)
        sys.exit(1)
    return file

#==== Data Locations ====
# Link data here or edit 
#========================
lc_epsg     = 'data/EPSG3413/EPSG3413grid.json'        #NOTE: created by plot_projections.py
lc_epsg_shr = 'data/EPSG3413/EPSG3413grid_shrunk.json' #NOTE: created by plot_grids.py
lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_seaRise  = 'data/SeaRise/Greenland1km.nc'
lc_racmo2p0 = 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc'
lc_InSAR    = 'data/InSAR/Joughin2015/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'
lc_mask     = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'


#==== SETUP ====
# get args, time
# load data sets
#===============
stamp = datetime.date.today().strftime("%Y_%m_%d")
f_base     = 'templates/greenland_1km.epsg3413.nc'
f_1km      = 'complete/greenland_1km_'+stamp+'.epsg3413.nc'
f_template = 'greenland.epsg3413.config'

# parse the command line arguments
parser = argparse.ArgumentParser()   # -h or --help automatically included!

parser.add_argument('-e', '--extended', help='Produce the extended grid.', action='store_true')

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

speak.notquiet(args,"\nBuilding the Greenland datasets in the EPSG:3413 projection.")
speak.notquiet(args,  "============================================================\n")

# load in datasets
speak.notquiet(args,"Loading the datasets.")

from data import epsg3413
f_epsg = abs_existing_file(lc_epsg)
speak.verbose(args,"   Found EPSG:3413 grid specs")
f_epsg_shr = abs_existing_file(lc_epsg_shr)
speak.verbose(args,"   Found shrunken EPSG:3413 grid specs")

from data import bamberdem
nc_bamber = get_nc_file(lc_bamber,'r')
speak.verbose(args,"   Found Bamber DEM")
    
from data import searise
nc_seaRise = get_nc_file(lc_seaRise,'r')
speak.verbose(args,"   Found Sea Rise data")

from data import racmo2p0
nc_racmo2p0 = get_nc_file(lc_racmo2p0,'r')
speak.verbose(args,"   Found RACMO 2.0 data")

from data import insar
try:
    nc_insar = get_nc_file(lc_InSAR,'r')
except Exception:
    speak.verbose(args,"\n   Building InSAR velocity dataset...\n")
    subprocess.call("python util/convert_velocities.py "+os.path.dirname(lc_InSAR), shell=True)
    nc_insar = get_nc_file(lc_InSAR,'r')
speak.verbose(args,"   Found InSAR data")

from data import icebridge
nc_massCon = get_nc_file(lc_massCon,'r')
speak.verbose(args,"   Found Mass Conserving Bed data")

from data import ice2sea
nc_mask = get_nc_file( lc_mask, 'r'  )
speak.verbose(args,"   Found Zurich mask")

speak.verbose(args,"\n   All data files found!")


#====== EPSG:3413 ======
# this is a 1km dataset 
#=======================
speak.notquiet(args,"\nCreating the base grid."),
nc_base, base = epsg3413.build_base(f_base, f_epsg, 1000.)
speak.notquiet(args,"   Done!")

#==== Projections ====
# All the projections 
# needed for the data 
#=====================
speak.notquiet(args,"\nGetting the projections.")
proj_epsg3413, proj_eigen_gl04c = projections.greenland(args, lc_bamber)
speak.notquiet(args,"   Done!")

#====== Lat,Lon ======
# Get the cell-center 
# lats and lons       
#=====================
speak.notquiet(args,"\nDetermining the  latitudes and longitudes of the grid cell-centers.")
projections.grid_center_latlons(nc_base, base, proj_epsg3413)
speak.notquiet(args,"   Done!")

#==== SeaRise Data =====
# this is a 1km dataset 
#=======================
speak.notquiet(args,"\nGetting bheatflx and presartm from the SeaRise data.")
searise.bheatflx_artm_epsg3413(args, nc_seaRise, nc_base, base, proj_epsg3413, proj_eigen_gl04c)
speak.notquiet(args,"   Done!")
nc_seaRise.close()

#==== RACMO2.0 Data =====
# this is a 1km dataset  
#========================
speak.notquiet(args,"\nGetting acab from the RACMO 2.0 data.")
racmo2p0.acab_epsg3413(args, nc_racmo2p0, nc_bamber, nc_base, base, proj_epsg3413, proj_eigen_gl04c)
speak.notquiet(args,"   Done!")
nc_racmo2p0.close()

#==== InSAR velocity Data ====
# this is a 500m dataset in   
# the ESPG-3413 projection    
#=============================
speak.notquiet(args,"\nGetting vy, vx, ey, and ex from the InSAR data.")
insar.velocity_epsg3413(args, nc_insar, nc_base, base)
speak.notquiet(args,"   Done!")
nc_insar.close()

#==== Mass Conserving Bed Data ===
# This is the new (2015) bed data 
#=================================
speak.notquiet(args,"\nGetting thk, topg, and topgerr from the mass conserving bed data.")
icebridge.mcb_epsg3413(args, nc_massCon, nc_bamber, nc_base, base, proj_epsg3413, proj_eigen_gl04c)
speak.notquiet(args,"   Done!")
nc_bamber.close()
nc_massCon.close()

#==== Done getting data ====
#===========================
nc_base.close()

#==== add time dim and shrink ====
# apply to all the variables and  
# shrink to size around ice sheet 
#=================================
speak.notquiet(args,"\nAdding the time dimension and creating the 1km dataset.")
epsg3413.add_time(args, f_base, f_1km, f_template, f_epsg_shr)

#==== Coarsen ==== 
# make 2, 4 and 8  
# km datasets      
#==================
speak.notquiet(args,"\nCreating coarser datasets.")
coarse_list = [2,4,5,8]   # in km
bamberdem.coarsen(args, f_1km, f_template, coarse_list)

#==== and done! ====
#===================
speak.notquiet(args,"\nFinished building the datasets.")
