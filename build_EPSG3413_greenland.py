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
lc_epsg     = 'EPSG3413grid.json'
lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_seaRise  = 'data/SeaRise/Greenland1km.nc'
lc_racmo2p0 = 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc'
#lc_InSAR    = 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_InSAR    = 'data/InSAR/Joughin2015/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'
lc_mask     = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'


#==== SETUP ====
# get args, time
# load data sets
#===============
stamp = datetime.date.today().strftime("%Y_%m_%d")
f_base = 'templates/greenland_1km.mcb.epsg3413.nc'

# parse the command line arguments
parser = argparse.ArgumentParser()   # -h or --help automatically included!

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
speak.verbose(args,"\nBuilding the base dataset: "+f_base)

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


#==== Mass Conserving Bed Data ===
# This is the new (2015) bed data 
#=================================
speak.notquiet(args,"\nGetting thk, topg, and topgerr from the mass conserving bed data.")

icebridge.mcb_epsg3413(args, nc_massCon, nc_base, base)

nc_bamber.close()
nc_massCon.close()
#==== Done getting data ====
#===========================
nc_base.close()
