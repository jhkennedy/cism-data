"""
data.racmo2p0 : RACMO 2.0 data import module.

This module provides functions to import data from the RACMO 2.0 
dataset into a CISM dataset. 

Functions list:
    * get_acab(args, nc_racmo2p0, nc_base, base)

Notes
-----
This data is from Ian Howat. No grid or projection information is
included with the NetCDF file, but it appears to use the same grid 
(albeit transposed) and projection as the Bamber DEM:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = 71 degrees
    * Latitude of projection origin = 90 degrees
    * Central meridian = -39 degrees
    * false eastings = 0
    * flase northings = 0
    * 1 km postings with
        + lower-left corner y,x: -3400000.0,-800000.0 (m) 
        + upper-right corner y,x:  600000.0, 700000.0 (m)

The data is on a 2501 x 3001 grid which, if transposed, is the same as the 
Bamber grid. These *should* be the same.
"""

import sys
import scipy
import pyproj
import numpy as np

from util.ncfunc import copy_atts, copy_atts_add_fill
from util import speak
from util import projections 

def acab_epsg3413(args, nc_racmo, nc_bamber, nc_base, base, proj_epsg3413, proj_eigen_gl04c):
    bamber = projections.DataGrid()
    bamber.y = nc_bamber.variables['projection_y_coordinate']
    bamber.x = nc_bamber.variables['projection_x_coordinate']
    bamber.ny = bamber.y[:].shape[0]
    bamber.nx = bamber.x[:].shape[0] 
    bamber.make_grid()
  

    speak.verbose(args,  '   Interpolating smb to acab in base grid.')
    sys.stdout.write("   [%-60s] %d%%" % ('='*0, 0.))
    sys.stdout.flush()
    racmo_data = nc_racmo.variables['smb'][:,::-1].transpose() / 910.
    racmo_data = np.ma.masked_invalid(racmo_data) # find invalid data and create a mask
    racmo_data = racmo_data.filled(0.)            # fill invalid data with zeros

    base2bamber = projections.transform(base, proj_epsg3413, proj_eigen_gl04c)
    racmo_interp = scipy.interpolate.RectBivariateSpline( bamber.y[:], bamber.x[:], racmo_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
    base_bamber = np.zeros( base.dims )
    for ii in range(0, base.nx):
        ctr = (ii*60)/base.nx
        sys.stdout.write("\r   [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
        sys.stdout.flush()
        base_bamber[:,ii] = racmo_interp.ev(base2bamber.y_grid[:,ii], base2bamber.x_grid[:,ii] )
    sys.stdout.write("\r   [%-60s] %d%%\n" % ('='*60, 100.))
    sys.stdout.flush()

    #NOTE: RectBivariateSpline extrapolates data outside the convex hull
    #      (constant value), so, we need to create a mask for values outisde
    #      the IceBridge convex hull...
    base_y_mask = np.ma.masked_outside(base2bamber.y_grid, bamber.y[0], bamber.y[-1])
    base_x_mask = np.ma.masked_outside(base2bamber.x_grid, bamber.x[0], bamber.x[-1])
    base_masked = np.ma.masked_array(base_bamber, mask=np.logical_or(base_y_mask.mask,base_x_mask.mask))
    
    base.var = nc_base.createVariable('acab', 'f4', ('y','x',) )
    base.var[:] = base_masked.filled(0.)  
    base.var.long_name = 'Water Equivalent Surface Mass Balance'
    base.var.standard_name = 'land_ice_lwe_surface_specific_mass_balance'
    base.var.units = 'mm/year' 
    base.var.grid_mapping = 'epsg_3413'
    base.var.coordinates = 'lon lat'
    base.var.source = 'Ian Howat' 
    base.var.comments = '1961--1990 mean surface mass balance from RACMO 2.0 '
    

def acab_bamber(args, nc_racmo2p0, nc_base, base):
    """Get acab from the RACMO 2.0 data.
 
    This function pulls in the `smb` variable from the RACMO 2.0 dataset 
    and writes it to the base dataset as `acab`. NetCDF attributes are preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    nc_racmo2p0 :
        An opened netCDF Dataset containing the RACMO 2.0 data.
    nc_base :
        The created netCDF Dataset that will contain the base data.
    base :
        A DataGrid() class instance that holds the base data grid information.
   
    """
    racmo2p0_data = np.ndarray( (base.ny, base.nx) )

    racmo2p0_data[:,:] = 0.
    racmo2p0_smb = nc_racmo2p0.variables['smb']
    racmo2p0_data[:,:] = racmo2p0_smb[:,::-1].transpose() / 910.
    racmo2p0_data = np.ma.masked_invalid(racmo2p0_data) # find invalid data and create a mask
    racmo2p0_data = racmo2p0_data.filled(0.)            # fill invalid data with zeros

    speak.verbose(args,"   Writing acab to base.")
    base_acab = nc_base.createVariable( 'acab', 'f4', ('y','x',) )
    base_acab[:,:] = racmo2p0_data[:,:]
    copy_atts(racmo2p0_smb, base_acab) #FIXME: check atribute units -- divided by 910 earlier

