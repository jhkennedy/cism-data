"""
data.ice2sea : Ice2Sea data import module.

This module provides functions to apply the Ice2Sea Greenland ice
mask, known as the Zurich mask, to a CISM dataset. 

Functions list:
    * apply_mask(args, nc_mask, nc_base)

Notes
-----
This data is associated with deliverable 3.2.8 of the Ice2Sea project.

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

The data is on the 3001 x 2501 Bamber DEM grid.

References
---------- 
Rastner, P., Bolch, T., Mölg, N., Machguth, H., Le Bris, R., and Paul, F.: The
first complete inventory of the local glaciers and ice caps on Greenland, The
Cryosphere, 6, 1483-1495, doi:10.5194/tc-6-1483-2012, 2012.
"""

from util import speak

def apply_mask(args, nc_mask, nc_base ):
    """Apply Zurich mask to thk and make usrf.

    This function pulls in the Zurich mask, applies it to the variable `thk`  
     and then creates the variable `usrf.`

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    nc_mask :
        An opened netCDF Dataset containing the Zurich mask.
    nc_base :
        The created netCDF Dataset that contains the base data.
    """
    base_thk  = nc_base.variables['thk']
    thk_data = base_thk[:,:]

    base_topg = nc_base.variables['topg']
    topg_data = base_topg[:,:]

    mask = nc_mask.variables['IceSheetMask']

    speak.verbose(args,"   Applying mask to thk.")
    thk_data = thk_data * mask[:,:]
    base_thk[:,:] = thk_data

    speak.verbose(args,"   Creating usrf.")
    base_usrf = nc_base.createVariable('usrf', 'f4',('y','x',))
    base_usrf[:,:] = thk_data + topg_data

