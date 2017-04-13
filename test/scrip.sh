#!/bin/bash
#
#
# Batch script to submit to create ESMF mapping file
#
# Set up for yellowstone
#
# yellowstone-specific batch commands:
#BSUB -P P93300601        # project number
#BSUB -n 8                # number of processors
#BSUB -R "span[ptile=16]"
#BSUB -W 1:00             # wall-clock limit
#BSUB -q caldera          # queue
#BSUB -o regrid.%J.out    # ouput filename
#BSUB -e regrid.%J.err    # error filename
#BSUB -J create_ESMF_map  # job name
#BSUB -N                  # send email upon job completion
 
#----------------------------------------------------------------------
 
#----------------------------------------------------------------------
# Set user-defined parameters here
#----------------------------------------------------------------------
 
# FIXME: Replace the following lines with paths to SCRIP grid files and names of your grids
# filesrc="/glade/p/cesmdata/cseg/inputdata/glc/cism/griddata/SCRIPgrid_greenland_4km_epsg3413_c161223.nc"
filesrc="/glade/p/cesmdata/cseg/inputdata/glc/cism/griddata/SCRIPgrid_greenland_4km_epsg3413_c170223.nc"
filedst="/glade/p/cesmdata/cseg/inputdata/glc/cism/griddata/SCRIPgrid_greenland_4km_epsg3413_c170223.nc"
namesrc='gland4kmOld'
namedst='gland4kmNew'
 
typesrc='regional'
typedst='regional'
maptype='aave'
 
#----------------------------------------------------------------------
# Done setting user-defined parameters
#----------------------------------------------------------------------
 
#----------------------------------------------------------------------
# Stuff done in a machine-specific way
#----------------------------------------------------------------------
 
# Determine number of processors we're running on
host_array=($LSB_HOSTS)
REGRID_PROC=${#host_array[@]}
 
#----------------------------------------------------------------------
# Begin general script
#----------------------------------------------------------------------
 
cmdargs="--filesrc $filesrc --filedst $filedst --namesrc $namesrc --namedst $namedst --typesrc $typesrc --typedst $typedst --maptype $maptype --batch"
env REGRID_PROC=$REGRID_PROC ./create_ESMF_map.sh $cmdargs
