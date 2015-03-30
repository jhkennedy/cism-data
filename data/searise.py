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

The data uses the same projection as the Bamber DEM:
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

import numpy as np

from util.ncfunc import copy_atts
from util import speak

def get_bheatflx_artm(args, nc_seaRise, nc_base, base) :
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

