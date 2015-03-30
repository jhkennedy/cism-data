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
Baber grid. These *should* be the same.
"""

import numpy as np

from util.ncfunc import copy_atts
from util import speak

def get_acab(args, nc_racmo2p0, nc_base, base):
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

