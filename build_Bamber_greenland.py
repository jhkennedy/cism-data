#!/usr/bin/env python

import os
import datetime
import subprocess
import argparse

from util import speak
from util import projections
from util.ncfunc import get_nc_file

"""
Build a CISM dataset
"""
#==== Data Locations ====
# Link data here or edit 
#========================
lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_seaRise  = 'data/SeaRise/Greenland1km.nc'
lc_racmo2p0 = 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc'
#lc_InSAR    = 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_InSAR    = 'data/InSAR/Joughin2015/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'
lc_mask     = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'


#==== SETUP =====
# get args, time 
# load data sets 
#================
stamp = datetime.date.today().strftime("%Y_%m_%d")
f_base = 'templates/greenland_1km.mcb.nc'

# parse the command line arguments
parser = argparse.ArgumentParser()   # -h or --help automatically included!

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

speak.notquiet(args,"\nBuilding the Greenland datasets in the Bamber projection.")
speak.notquiet(args,  "=========================================================\n")

# load in datasets
speak.notquiet(args,"Loading the datasets.")

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


#===== Bamber DEM ======
# this is a 1km dataset 
#=======================
speak.verbose(args,"\nBuilding the base dataset: "+f_base)

speak.notquiet(args,"\nCreating the base grid."),

nc_base, base = bamberdem.build_base(f_base, nc_bamber)

speak.notquiet(args,"   Done!")

#==== Projections ====
# All the projections 
# needed for the data 
#=====================
speak.notquiet(args,"\nGetting the projections.")

proj_epsg3413, proj_eigen_gl04c = projections.greenland(args, lc_bamber)

speak.notquiet(args,"   Done!")

# transform meshes. 
speak.verbose(args,"   Creating the transform meshes: base Bamber grid to EPSG-3413.")

trans = projections.transform(base, proj_eigen_gl04c, proj_epsg3413)

speak.notquiet(args,"   Done!")

#==== SeaRise Data =====
# this is a 1km dataset 
#=======================
speak.notquiet(args,"\nGetting bheatflx and presartm from the SeaRise data.")

searise.get_bheatflx_artm(args, nc_seaRise, nc_base, base)

nc_seaRise.close()
#==== RACMO2.0 Data =====
# this is a 1km dataset  
#========================
speak.notquiet(args,"\nGetting acab from the RACMO 2.0 data.")

racmo2p0.get_acab(args, nc_racmo2p0, nc_base, base)

nc_racmo2p0.close()
#==== InSAR velocity Data ====
# this is a 500m dataset in   
# the ESPG-3413 projection    
#=============================
speak.notquiet(args,"\nGetting vy, vx, ey, and ex from the InSAR data.")

insar.get_velocity(args, nc_insar, nc_base, trans)

nc_insar.close()
#==== Mass Conserving Bed Data ===
# This is the new (2015) bed data 
#=================================
speak.notquiet(args,"\nGetting thk, topg, and topgerr from the mass conserving bed data.")

icebridge.mcb_bamber(args, nc_massCon, nc_bamber, nc_base, base, trans, proj_eigen_gl04c, proj_epsg3413)

nc_bamber.close()
nc_massCon.close()
#==== Zurich mask =====
# apply mask, and get  
# new surface variable 
#======================
speak.notquiet(args,"\nGetting the Zurich Mask.")

base = None
nc_base.close()   # need to read in some data from nc_base now
nc_base = get_nc_file(f_base,'r+')

ice2sea.apply_mask(args, nc_mask, nc_base)

nc_mask.close()
#==== Done getting data ====
#===========================
nc_base.close()

#==== add time dim and shrink ====
# apply to all the variables and  
# shrink to size around ice sheet 
#=================================
speak.notquiet(args,"\nAdding the time dimension and creating the 1km dataset.")

f_1km      = 'complete/greenland_1km_'+stamp+'.mcb.nc'
f_template = 'greenland.mcb.config'

bamberdem.add_time(args, f_base, f_1km, f_template)

#==== Coarsen ==== 
# make 2, 4 and 8  
# km datasets      
#==================
speak.notquiet(args,"\nCreating coarser datasets.")

coarse_list = [2,4,8]   # in km

bamberdem.coarsen(args, f_1km, f_template, coarse_list)

#==== and done! ====
#===================
speak.notquiet(args,"\nFinished building the datasets.")
