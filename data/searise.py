"""
data.searise : Sea Rise data import module.

This module provides functions to import data from the Sea Rise dataset
into a CISM dataset. 

Functions list:
    * get_bheatflx_artm(args, nc_seaRise, nc_base, base)

Notes
-----
This data is associated with the searise project. More information
can be found at:

http://websrv.cs.umt.edu/isis/index.php/1km_Greenland_data_set

The data uses the Bamber DEM projection:
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

The data is on a 2801 x 1501 grid, with the corresponding lower-left
bamber point being [100,500].

References
----------
Bamber, J.L., R.L. Layberry, S.P. Gogenini. 2001. A new ice thickness 
and bed data set for the Greenland ice sheet 1: Measurement, data 
reduction, and errors. Journal of Geophysical Research 106 (D24): 
33773-33780. 

Bamber, J.L., R.L. Layberry, S.P. Gogenini. 2001. A new ice thickness 
and bed data set for the Greenland ice sheet 2: Relationship between 
dynamics and basal topography. Journal of Geophysical Research 106 (D24): 
33781-33788.
"""

import sys
import scipy
import pyproj
import numpy as np

from util.ncfunc import copy_atts, copy_atts_add_fill
from util import speak
from util import projections 

def bheatflx_artm_epsg3413(args, nc_seaRise, nc_base, base, proj_epsg3413, proj_eigen_gl04c):
    seaRise = projections.DataGrid()
    seaRise.y = nc_seaRise.variables['y']
    seaRise.x = nc_seaRise.variables['x']
    seaRise.ny = seaRise.y[:].shape[0]
    seaRise.nx = seaRise.x[:].shape[0]
    seaRise.make_grid()

    base_vars = {'bheatflx':'bheatflx', 'artm':'presartm'}
    for base_var, sea_var in base_vars.iteritems():
        speak.verbose(args,'   Interpolating '+sea_var+' to '+base_var+' in base grid.')
        sys.stdout.write("   [%-60s] %d%%" % ('='*0, 0.))
        sys.stdout.flush()
        sea_data = nc_seaRise.variables[sea_var][0,:,:]

        base2bamber = projections.transform(base, proj_epsg3413, proj_eigen_gl04c)
        seaRise_interp = scipy.interpolate.RectBivariateSpline( seaRise.y[:], seaRise.x[:], sea_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
        base_bamber = np.zeros( base.dims )
        for ii in range(0, base.nx):
            ctr = (ii*60)/base.nx
            sys.stdout.write("\r   [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
            sys.stdout.flush()
            base_bamber[:,ii] = seaRise_interp.ev(base2bamber.y_grid[:,ii], base2bamber.x_grid[:,ii] )
        sys.stdout.write("\r   [%-60s] %d%%\n" % ('='*60, 100.))
        sys.stdout.flush()

        #NOTE: RectBivariateSpline extrapolates data outside the convex hull
        #      (constant value), so, we need to create a mask for values outisde
        #      the IceBridge convex hull...
        base_y_mask = np.ma.masked_outside(base2bamber.y_grid, seaRise.y[0], seaRise.y[-1])
        base_x_mask = np.ma.masked_outside(base2bamber.x_grid, seaRise.x[0], seaRise.x[-1])
        base_masked = np.ma.masked_array(base_bamber, mask=np.logical_or(base_y_mask.mask,base_x_mask.mask))

        if base_var != 'bheatflx':
            base_bamber[base_masked.mask] = -9999.

        base.var = nc_base.createVariable(base_var, 'f4', ('y','x',) )
        if base_var == 'bheatflx':
            base.var[:] = -base_bamber[:] # invert sign!
        else:
            base.var[:] = base_bamber[:]

        copy_atts_add_fill(nc_seaRise.variables[sea_var], base.var, -9999.)
        base.var.grid_mapping = 'epsg_3413'
        base.var.coordinates = 'lon lat'
     
def bheatflx_artm_bamber(args, nc_seaRise, nc_base, base) :
    """Get bheatflx and artm from the sea rise data.

    This function pulls in the `bheatflx` and `presartm` variables from the
    Sea Rise dataset and writes them to the base dataset as `bheatflx` and
    `artm`. NetCDF attributes are preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    nc_seaRise :
        An opened netCDF Dataset containing the Sea Rise data.
    nc_base :
        The created netCDF Dataset that will contain the base data.
    base :
        A DataGrid() class instance that holds the base data grid information.
    """
    seaRise_y = nc_seaRise.variables['y']
    seaRise_ny = seaRise_y[:].shape[0]

    seaRise_x = nc_seaRise.variables['x']
    seaRise_nx = seaRise_x[:].shape[0]

    # to convert to Bamber 1km grid
    seaRise_data = np.ndarray( (base.ny, base.nx) )
    seaRise_y_equal_base = 100
    seaRise_x_equal_base = 500

    # get basal heat flux
    #--------------------
    seaRise_data[:,:] = 0.
    seaRise_bheatflx = nc_seaRise.variables['bheatflx']
    seaRise_data[ seaRise_y_equal_base:seaRise_y_equal_base+seaRise_ny , seaRise_x_equal_base:seaRise_x_equal_base+seaRise_nx ] = -seaRise_bheatflx[0,:,:] # invert sign!

    speak.verbose(args,"   Writing bheatflx to base.")
    base_bheatflx = nc_base.createVariable('bheatflx', 'f4', ('y','x',) )
    base_bheatflx[:,:] = seaRise_data[:,:] 
    copy_atts(seaRise_bheatflx, base_bheatflx)

    # get annual mean air temperature (2m)
    #-------------------------------------
    seaRise_data[:,:] = 0.
    seaRise_presartm = nc_seaRise.variables['presartm']
    seaRise_data[ seaRise_y_equal_base:seaRise_y_equal_base+seaRise_ny , seaRise_x_equal_base:seaRise_x_equal_base+seaRise_nx ] = seaRise_presartm[0,:,:]

    speak.verbose(args,"   Writing artm to base.")
    base_artm = nc_base.createVariable('artm', 'f4', ('y','x',) )
    base_artm[:,:] = seaRise_data[:,:] 
    copy_atts(seaRise_presartm, base_artm)

