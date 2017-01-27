#!/usr/bin/env python

import os
import json
import argparse
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt

from util import speak
from util import projections
from util.ncfunc import get_nc_file


#lc_base   = 'complete/greenland_1km_2016_11_23.mcb.nc'
lc_base   = 'complete/greenland_1km_2016_03_17.mcb.nc'


#==================================
# parse the command line arguments 
#==================================
parser = argparse.ArgumentParser()   # -h or --help automatically included!

volume = parser.add_mutually_exclusive_group()
volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

args = parser.parse_args()

speak.notquiet(args,"\nPlotting the representation of the Bamber grid in the EPSG:3413 projection")
speak.notquiet(args,  "==========================================================================\n")


#==================
# load in datasets 
#==================
speak.notquiet(args,"Loading the datasets.")

nc_base = get_nc_file(lc_base,'r')
speak.verbose(args,"   Found Bamber DEM")

speak.verbose(args,"\n   All data files found!")


#===== Bamber DEM ======
# this is a 1km dataset 
#=======================

speak.notquiet(args,"\nCreating the 1 km Bamber grid."),
base = projections.DataGrid()
base.y = nc_base.variables['y1']
base.x = nc_base.variables['x1']
base.ny = base.y[:].shape[0]
base.nx = base.x[:].shape[0] 
base.make_grid()

speak.notquiet(args,"   Done!")


#==== Projections ====
# All the projections 
# needed for the data 
#=====================
speak.notquiet(args,"\nGetting the projections.")

proj_epsg3413, proj_eigen_gl04c = projections.greenland()

speak.notquiet(args,"   Done!")


#========================
# Determine the new grid 
#========================
speak.notquiet(args, "\nProject the Bamber grid into EPSG:3413.")

B_ll = ( np.amin(base.x_grid[1,:]), np.amin(base.y_grid[1,:]) )
B_lr = ( np.amax(base.x_grid[1,:]), np.amax(base.y_grid[1,:]) )
B_ur = ( np.amax(base.x_grid[-1,:]), np.amax(base.y_grid[-1,:]) )
B_ul = ( np.amin(base.x_grid[-1,:]), np.amin(base.y_grid[-1,:]) )

B_xs = np.array([ B_ll[0], B_lr[0], B_ur[0], B_ul[0], B_ll[0] ])
B_ys = np.array([ B_ll[1], B_lr[1], B_ur[1], B_ul[1], B_ll[1] ])

B_lons, B_lats = proj_eigen_gl04c(B_xs, B_ys, inverse=True) 

B2E_xs, B2E_ys = proj_epsg3413(B_lons, B_lats)

speak.verbose(args,  "    Bamber in EPSG:3413:")
speak.verbose(args,  "      Lower Left  (x,y): ("+str(B2E_xs[0])+", "+str(B2E_ys[0])+")")
speak.verbose(args,"\n      Lower Right (x,y): ("+str(B2E_xs[1])+", "+str(B2E_ys[1])+")")
speak.verbose(args,"\n      Upper Right (x,y): ("+str(B2E_xs[2])+", "+str(B2E_ys[2])+")")
speak.verbose(args,"\n      Upper Left  (x,y): ("+str(B2E_xs[3])+", "+str(B2E_ys[3])+")")

speak.notquiet(args, "\nDeterniming the new EPSG:3413 grid from the transformed Bamber grid.")

e_ll = [np.floor(np.mean([B2E_xs[0], B2E_xs[3]])/1000.)*1000., 
        np.floor(np.mean([B2E_ys[0], B2E_ys[1]])/1000.)*1000.]

e_ur = [ np.ceil(np.mean([B2E_xs[1], B2E_xs[2]])/1000.)*1000., 
         np.ceil(np.mean([B2E_ys[3], B2E_ys[2]])/1000.)*1000.]


E_xs = np.array([e_ll[0], e_ur[0], e_ur[0], e_ll[0], e_ll[0]])
E_ys = np.array([e_ll[1], e_ll[1], e_ur[1], e_ur[1], e_ll[1]])

E_lons, E_lats = proj_epsg3413(E_xs, E_ys, inverse=True)

offset = -80000.0

O_xs = np.array([e_ll[0]+offset, e_ur[0], e_ur[0], e_ll[0]+offset, e_ll[0]+offset])
O_ys = np.array([e_ll[1], e_ll[1], e_ur[1], e_ur[1], e_ll[1]])

O_lons, O_lats = proj_epsg3413(O_xs, O_ys, inverse=True)

speak.notquiet(args,  "    New EPSG:3413 grid:")
speak.notquiet(args,  "      Lower Left  (x,y): ("+str(E_xs[0])+", "+str(E_ys[0])+")")
speak.notquiet(args,"\n      Lower Right (x,y): ("+str(E_xs[1])+", "+str(E_ys[1])+")")
speak.notquiet(args,"\n      Upper Right (x,y): ("+str(E_xs[2])+", "+str(E_ys[2])+")")
speak.notquiet(args,"\n      Upper Left  (x,y): ("+str(E_xs[3])+", "+str(E_ys[3])+")")


#=====================================
# Write the grid specs to a json file 
#=====================================
O_grid = {'ll':[O_xs[0], O_ys[0]],
          'ur':[O_xs[2], O_ys[2]],
          'xs':O_xs.tolist(),
          'ys':O_ys.tolist(),
          'lons': O_lons.tolist(),
          'lats': O_lats.tolist(),
          'projstring':'+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m'
         }

with open('EPSG3413grid_shrunk.json', 'w') as f:
    json.dump(O_grid, f)


#======
# Plot 
#======
speak.notquiet(args,"\nPlotting the grid.")

plt.figure(1, figsize=(10,8), dpi=150)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


#======= Plot ========
# EPSG:3413 in Bamber 
#=====================
plt.subplot(1,2,1)

#NOTE: Basemap adds false Eastings and Northings dependent upon your corners.
#      The proj4 projections, however, don't have any, so you'll alsways have
#      to passs the map longs and lats. 
map_ll_lon, map_ll_lat = proj_eigen_gl04c(B_xs[0]-300000., B_ys[0]-600000., inverse=True) 
map_ur_lon, map_ur_lat = proj_eigen_gl04c(B_xs[2]+300000., B_ys[2]+300000., inverse=True) 
glmap = Basemap(llcrnrlon=map_ll_lon, llcrnrlat=map_ll_lat, 
                urcrnrlon=map_ur_lon, urcrnrlat=map_ur_lat,
                projection='stere', lat_ts=71.0, lat_0=90.0, lon_0=321.0, 
                resolution='l')

glmap.fillcontinents(color='gray', lake_color='white')
glmap.drawcoastlines()

B_X, B_Y = glmap(B_lons, B_lats)
glmap.plot(B_X, B_Y, 'bo-', label='Bamber grid')

E2B_X, E2B_Y = glmap(E_lons, E_lats)
glmap.plot(E2B_X, E2B_Y, 'mo-', label='EPSG:3413 grid')

O2B_X, O2B_Y = glmap(O_lons, O_lats)
glmap.plot(O2B_X, O2B_Y, 'go-', label='Offset EPSG:3413 grid')

plt.title('Greenland using the Bamber projection')
plt.legend(loc='lower center', ncol=2, fancybox=True, shadow=True)


#======= Plot ========
# Bamber in EPSG:3413 
#=====================
plt.subplot(1,2,2)

#NOTE: Basemap adds false Eastings and Northings dependent upon your corners.
#      The proj4 projections, however, don't have any, so you'll alsways have
#      to passs the map longs and lats. 
map_ll_lon, map_ll_lat = proj_epsg3413(E_xs[0]-300000., E_ys[0]-600000., inverse=True) 
map_ur_lon, map_ur_lat = proj_epsg3413(E_xs[2]+300000., E_ys[2]+300000., inverse=True) 
glmap = Basemap(llcrnrlon=map_ll_lon, llcrnrlat=map_ll_lat,
                urcrnrlon=map_ur_lon, urcrnrlat=map_ur_lat,
                resolution='l', epsg=3413)

glmap.fillcontinents(color='gray', lake_color='white')
glmap.drawcoastlines()

B2E_X, B2E_Y = glmap(B_lons, B_lats)
glmap.plot(B2E_X, B2E_Y, 'bo-', label='Bamber grid')

E_X, E_Y = glmap(E_lons, E_lats)
glmap.plot(E_X, E_Y, 'mo-', label='EPSG:3413 grid')

O_X, O_Y = glmap(O_lons, O_lats)
glmap.plot(O_X, O_Y, 'go-', label='Offset EPSG:3413 grid')

plt.title('Greenland using the EPSG:3413 projection')
plt.legend(loc='lower center', ncol=2, fancybox=True, shadow=True)


#========
# Finish 
#========
plt.tight_layout()
plt.show()

nc_base.close()
speak.notquiet(args,"   Done!\n")



