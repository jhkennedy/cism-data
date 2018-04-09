#!/bin/csh

# load ncl to have weight-generation command ESMF_RegridWeightGen used in ncremap
module load ncl/6.1.0
module load nco


set cesmfile  = /lustre/atlas1/cli106/proj-shared/4ue/CESM-CISM/b.e10.BG20TRCN.f09_g16.002/atm/climos/b.e10.BG20TRCN.f09_g16.002_DJF_195912_200502_climo.nc
set racmofile = /lustre/atlas1/cli106/proj-shared/4ue/CESM-CISM/racmo23_GRN_monthly/climos_1960-2005/racmo23_GRN_monthly.t2m.1960-2005.DJF.nc 
set racmolatlon = /lustre/atlas1/cli106/proj-shared/4ue/CESM-CISM/racmo23_GRN_monthly/RACMO23_masks_ZGRN11.nc
set tmpfile1 = racmo23_1960-2005.DJF_latlont2m.nc 
set tmpfile2 = cesm_DJF_latlonTREFHT.nc 
set ReMapFile = ReMap_b.e10.BG20TRCN.f09_g16.002_DJF_195912_200502_climo.nc 

##-- remap cesm to have same size as racmo
##-first, append lat and lon from racmolatlon file to racmofile
ncks -C -v lat,lon $racmolatlon -O $tmpfile1
ncks -C -v t2m $racmofile -A $tmpfile1

##-save the variables TREFHT from the cesmfile to a new file for remapping
ncks -C -v lat,lon,TREFHT $cesmfile -O $tmpfile2

##-- remap the variable TREFHT in cesm to have the same grid size as variable t2m in racmo
ncremap -v TREFHT -V t2m -a bilinear -i $tmpfile2 -d $tmpfile1 -m map.nc -o $ReMapFile










