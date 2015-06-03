"""
data.insar : InSAR data import module.

This module provides functions to import data from the MEaSUREs Greenland ice 
sheet velocity map from InSAR dataset into a CISM dataset. 

Functions list:
    * get_velocity(args, nc_insar, nc_base, trans)

Notes
-----
This data set, part of the NASA Making Earth System Data Records for Use in 
Research Environments (MEaSUREs) program, provides annual ice-sheet-wide 
velocity maps for Greenland, derived using Interferometric Synthetic Aperture 
Radar (InSAR) data from the RADARSAT-1 satellite. The data set currently 
contains ice velocity data for the winter of 2000-2001 and 2005-2006, 2006-2007,
and 2007-2008 acquired from RADARSAT-1 InSAR data from the Alaska Satellite 
Facility (ASF), and a 2008-2009 mosaic derived from the Advanced Land Observation 
Satelitte (ALOS) and TerraSAR-X data.

More information can be found at:

http://nsidc.org/data/docs/measures/nsidc0478_joughin/index.html

The data uses the ESPG:3413 projection which is different than the Bamber DEM:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = 70 degrees
    * Latitude of projection origin = 90 degrees
    * Central meridian = -45 degrees
    * false eastings = 0
    * flase northings = 0
    * 500 m postings with
        + lower-left corner y,x: -3370000.0,-645000.0 (m) 
        + upper-right corner y,x: -640500.0, 859500.0 (m)

The data is provided as a set of geodat files which are converted into a NetCDF 
dataset with variables vx, vy, ex, and ey. The NetCDF dataset is on a 5460 x 
3010 grid.

References
----------
Joughin, I., B. Smith, I. Howat, and T. Scambos.. 2010. MEaSUREs Greenland Ice 
Sheet Velocity Map from InSAR Data. [indicate subset used]. Boulder, Colorado 
USA: National Snow and Ice Data Center. 

Papers discussing this data include:

Joughin, I., B. Smith, I. Howat, T. Scambos, and T. Moon. 2010. Greenland Flow 
Variability from Ice-Sheet-Wide Velocity Mapping. Journal of Glaciology 
56(197): 415-430. 

Moon, T., Joughin, I., Smith, B. & Howat, I. 21st-Century Evolution of Greenland 
Outlet Glacier Velocities. Science 336, 576-578 (2012).

2012 Set
--------
As described above.

2015 Set
--------
The 2015 dataset contains an updated InSAR velocity map from I. Joughin (sent 
April, 2015). This is PRELIMINARY and UNPUBLISHED data. It attempts to fill some
holes and correct/reduce errors in the interior. 

"""

import numpy as np
from scipy import interpolate

from util.ncfunc import copy_atts
from util import speak

def get_velocity(args, nc_insar, nc_base, trans):
    """Get the velocities from the insar data.

    This function pulls in the `vx`, `vy`, `ex`  and `ey` variables from the
    InSAR dataset and writes them to the base dataset. NetCDF attributes are 
    preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    nc_insar :
        An opened netCDF Dataset containing the InSAR data.
    nc_base :
        The created netCDF Dataset that will contain the base data.
    trans :
        A DataGrid() class instance that holds the base data grid transformed 
        to the EPSG:3413 projection.
    """
    insar_y = nc_insar.variables['y']
    insar_ny = insar_y[:].shape[0]

    insar_x = nc_insar.variables['x']
    insar_nx = insar_x[:].shape[0]

    insar_data = np.ndarray( (insar_ny, insar_nx) )
    base_data = np.ndarray( (trans.ny,trans.nx) )


    for vv in ['vy','vx','ey','ex'] :
        insar_data[:,:] = 0.
        base_data[:,:] = 0.
        
        insar_var = nc_insar.variables[ vv ]
        insar_data[:,:] = insar_var[:,:]

        speak.verbose(args,"   Interpolating "+vv+".")
        insar_to_base = interpolate.RectBivariateSpline( insar_y[:], insar_x[:], insar_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster

        for ii in range(0, trans.nx):
            base_data[:,ii] = insar_to_base.ev(trans.y_grid[:,ii], trans.x_grid[:,ii] )
        
        data_min = np.amin( insar_data[ insar_data[:,:] != -2.e9] )
        data_max = np.amax( insar_data[ insar_data[:,:] != -2.e9] )

        data_mask = np.ma.masked_outside(base_data, data_min, data_max)
        base_data[ data_mask.mask ] = -2.e9
        
        speak.verbose(args,"   Writing "+vv+" to base.")
        base_var = nc_base.createVariable( vv, 'f4', ('y','x',) )
        base_var[:,:] = base_data[:,:]
        copy_atts(insar_var, base_var)

