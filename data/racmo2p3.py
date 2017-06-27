"""
data.racmo2p3 : RACMO 2.3 data import module.

This module provides functions to import data from a 1km downscaled RACMO 2.3 
dataset into a CISM dataset. 

Functions list:
    * acab_epsg3413(args, nc_racmo, nc_base, base)

Notes
-----
This data was provided by Jan Lenaerts from the citation below.  

The data uses the ESPG:3413 projection:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = 70 degrees
    * Latitude of projection origin = 90 degrees
    * Central meridian = -45 degrees
    * false eastings = 0
    * flase northings = 0


References
----------
Noel, B., van de Berg, W. J., Machguth, H., Lhermitte, S., Howat, I., Fettweis,
X., and van den Broeke, M. R.: a daily, 1 km resolution data set of downscaled
Greenland ice sheet surface mass balance (1958--2015), The Cryosphere, 10,
2361-2377, doi:10.5194/tc-10-2361-2016, 2016. 
"""

import sys
import scipy
import numpy as np

from util.ncfunc import copy_atts_bad_fill
from util import speak
from util import projections 

def acab_epsg3413(args, nc_racmo, nc_base, base):
    racmo = projections.DataGrid()
    racmo.y = nc_racmo.variables['y']
    racmo.x = nc_racmo.variables['x']
    racmo.ny = racmo.y[:].shape[0]
    racmo.nx = racmo.x[:].shape[0]
    racmo.make_grid()
    
    base_vars = {'smb':'SMB_rec'}
    for bvar, rvar in base_vars.iteritems():
        speak.verbose(args,'   Interpolating '+bvar+' and writing to base.')
        sys.stdout.write("   [%-60s] %d%%" % ('='*0, 0.))
        sys.stdout.flush()
        racmo_data = np.ma.masked_values(nc_racmo.variables[rvar][-1,:,:], 9.96921e+36)
        data_min = racmo_data.min() 
        data_max = racmo_data.max() 

        racmo_to_base = scipy.interpolate.RectBivariateSpline( racmo.y[::-1], racmo.x[:], racmo_data[::-1,:], kx=1, ky=1, s=0) # regular 2d linear interp. but faster
        base_data = np.zeros( base.dims )
        for ii in range(0, base.nx):
            ctr = (ii*60)/base.nx
            sys.stdout.write("\r   [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
            sys.stdout.flush()
            base_data[:,ii] = racmo_to_base.ev(base.y_grid[:,ii], base.x_grid[:,ii] )
        sys.stdout.write("\r   [%-60s] %d%%\n" % ('='*60, 100.))
        sys.stdout.flush()
        
        base_data[base_data < data_min] = 9.96921e+36
        base_data[base_data > data_max] = 9.96921e+36
        
        base.var = nc_base.createVariable(bvar, 'f4', ('y','x',) )
        base.var[:] = base_data[:]  
        copy_atts_bad_fill(nc_racmo.variables[rvar], base.var, 9.96921e+36)
        base.var.long_name = 'Water Equivalent Surface Mass Balance'
        base.var.standard_name = 'land_ice_lwe_surface_specific_mass_balance'
        base.var.units = 'mm/year'
        base.var.grid_mapping = 'epsg_3413'
        base.var.coordinates = 'lon lat'
        base.var.source = 'Jan Lenaerts'
        base.var.reference = 'Noel, B., van de Berg, W. J., Machguth, H., Lhermitte, S., Howat, I., Fettweis, X., and van den Broeke, M. R.: A daily, 1 km resolution data set of downscaled Greenland ice sheet surface mass balance (1958--2015), The Cryosphere, 10, 2361-2377, doi:10.5194/tc-10-2361-2016, 2016.'




