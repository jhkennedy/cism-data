"""
data.bamber : bamberdem data import module.

This module provides functions to import data from the Bamber dataset
and build a base dataset in the Bamber DEM. 

Functions list:
    * build_base (f_base, nc_bamber):
    * add_time (args, f_base, f_1km, f_template):
    * coarsen(args, f_base, f_template, coarse_list):

Notes
-----
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

The data is on a 3001 x 2501 grid.

References
----------
Bamber, J.L. et al., A new bed elevation dataset for Greenland, The Cryosphere, 
7, 499-510, doi:10.5194/tc-7-499-2013, 2013.
"""

from netCDF4 import Dataset

from util.ncfunc import copy_atts
from util.projections import DataGrid


def build_base (f_base, nc_bamber):
    """Build the Bamber base grid.

    This function opens the Bamber DEM dataset, pulls in the 
    `projection_y_coordinate` and `projection_x_coordinate` variables, creates 
    a new base dataset, creates a DataGrid() class instance to hold the base 
    grid, and sets up the base coordinate dimension and variables `y` and `x` 
    in the new base dataset. NetCDF attributes are preserved.

    Parameters
    ----------
    f_base :
        Filename for the created netCDF Dataset that will contain the base data.
    nc_bamber :
        An opened netCDF Dataset containing the Bamber dataset.
    """
    nc_base = Dataset(f_base,'w', format='NETCDF4')

    bamber_y = nc_bamber.variables['projection_y_coordinate']
    bamber_ny = bamber_y[:].shape[0] # number of y points for 1km grid

    bamber_x = nc_bamber.variables['projection_x_coordinate']
    bamber_nx = bamber_x[:].shape[0] # number of x points for 1km grid

    # make bamber 1km grid for base
    base = DataGrid()
    base.ny = bamber_ny
    base.nx = bamber_nx

    # create our base dimensions
    nc_base.createDimension('y',base.ny)
    nc_base.createDimension('x',base.nx)

    # create some base variables
    base.y = nc_base.createVariable('y', 'f4', 'y')
    base.y[:] = bamber_y[:]
    copy_atts(bamber_y, base.y)  # FIXME: units say km, but it's actuall in m

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = bamber_x[:]
    copy_atts(bamber_x, base.x)  # FIXME: units say km, but it's actuall in m

    # create some grids for interpolation
    base.make_grid()

    return (nc_base, base)
