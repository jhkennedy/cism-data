#!/usr/bin/env python2

import os
import sys
import datetime
import argparse

from util import speak
from util import finalize

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

#==== SETUP ====
# get args, time
# load data sets
#===============
stamp = datetime.date.today().strftime("%Y_%m_%d")
f_base     = 'templates/greenland_1km.epsg3413.nc'
f_1km      = 'complete/greenland_1km_'+stamp+'.epsg3413.nc'
f_template = 'templates/greenland.epsg3413.config'

from data import epsg3413
lc_epsg_shr = 'data/EPSG3413/EPSG3413grid_shrunk.json' #NOTE: created by plot_grids.py
f_epsg_shr = abs_existing_file(lc_epsg_shr)

from data import bamberdem

#==== add time dim and shrink ====
# apply to all the variables and  
# shrink to size around ice sheet 
#=================================
speak.notquiet(args,"\nAdding the time dimension and creating the 1km dataset.")
finalize.add_time_and_shrink(args, 'epsg_3413', f_base, f_1km, f_template, f_epsg_shr)

#==== Coarsen ==== 
# make 2, 4 and 8  
# km datasets      
#==================
speak.notquiet(args,"\nCreating coarser datasets.")
coarse_list = [2,4,5,8]   # in km
finalize.coarsen(args, 'epsg_3413', f_1km, f_template, coarse_list)


