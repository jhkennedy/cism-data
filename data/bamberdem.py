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
The data uses the same projection:
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

import os
from netCDF4 import Dataset

from util import speak
from util.ncfunc import copy_atts, get_nc_file
from util.projections import DataGrid
from templates import config

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
    nc_base = Dataset(f_base,'w')

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
    copy_atts(bamber_y, base.y) #FIXME: units say km, but it's actuall in m

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = bamber_x[:]
    copy_atts(bamber_x, base.x) #FIXME: units say km, but it's actuall in m

    # create some grids for interpolation
    base.make_grid()

    return (nc_base, base)

def add_time (args, f_base, f_1km, f_template):
    """Add the time dimension to a Bamber 1km DEM dataset and write config files. 

    This function opens the base dataset, creates a new 1km dataset with the time
    dimension that has been shrunken to just around the ice sheet, and creates a 
    CISM config file. NetCDF attributes are preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    f_base :
        Filename for the created netCDF Dataset that contains the base data.
    f_1km :
        Filename for the created netCDF Dataset that will contain the 1 km dataset.
    f_template :
        Filename for the template used to create the CISM config file.
    """
    # shrink dataset to the ice sheet
    nc_base = Dataset(f_base,'r')
    
    y_shrink = [100,2900+1] #NOTE: python stop exclusive, nco stop inclusive!
    x_shrink = [500,2000+1] #NOTE: python stop exclusive, nco stop inclusive!

    base = DataGrid()
    base.y = nc_base.variables['y']
    base.x = nc_base.variables['x']
    if args.extended:
        base.ny = base.y.shape[0]
        base.nx = base.x.shape[0]
    else:
        base.ny = base.y[ y_shrink[0]:y_shrink[1] ].shape[0]
        base.nx = base.x[ x_shrink[0]:x_shrink[1] ].shape[0]

    speak.verbose(args,"   Writing "+f_1km)
    nc_1km = Dataset(f_1km, 'w')
    nc_1km.createDimension('time', None )
    nc_1km.createDimension('y1', base.ny)
    nc_1km.createDimension('x1', base.nx)

    time = nc_1km.createVariable('time', 'f4', ('time',))
    y1   = nc_1km.createVariable('y1',   'f4', ('y1',)  )
    x1   = nc_1km.createVariable('x1',   'f4', ('x1',)  )

    copy_atts(base.y, y1)
    copy_atts(base.x, x1)

    if args.extended:
        y1[:] = base.y[:]
        x1[:] = base.x[:]
    else:
        y1[:] = base.y[ y_shrink[0]:y_shrink[1] ]
        x1[:] = base.x[ x_shrink[0]:x_shrink[1] ]
    time[0] = 0.

    for var_name, var_data in nc_base.variables.iteritems() : 
        if (var_name != 'x' and var_name != 'y'):
            var_1km = nc_1km.createVariable(var_name, 'f4', ('time','y1','x1',))
            if args.extended:
                var_1km[0,:,:] = var_data[:,:]
            else:
                var_1km[0,:,:] = var_data[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]
            copy_atts(var_data, var_1km)
    #copy_atts(var_data, var_1km)

    nc_base.close()
    nc_1km.close()
    os.chmod(f_1km, 0o644)   # uses an octal number!

    speak.verbose(args,"   Writing the 1km config file.")

    config.write(f_1km, f_template, base, 1)


def coarsen(args, f_base, f_template, coarse_list):
    """Coarsen an base dataset.

    This function opens the 1 km dataset and creates coarser resolution datasets.
    All NetCDF data attributes are preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    f_base :
        Filename for the base dataset to coarsen.
    f_template :
        Filename for the template used to create the CISM config files.
    coarse_list :
        List of resolutions to coarsen to. 
    """
    nc_base = Dataset(f_base,'r' )
    
    base = DataGrid()
     
    base.y = nc_base.variables['y1']
    base.ny = base.y[:].shape[0]

    base.x = nc_base.variables['x1']
    base.nx = base.x[:].shape[0]

    coarse_names = []
    for skip in coarse_list :
        coarse = DataGrid()

        coarse.ny = base.y[::skip].shape[0]
        coarse.nx = base.x[::skip].shape[0]
       
        idx = f_base.find('km')
        f_coarse = f_base[:idx-1]+str(skip)+f_base[idx:]

        speak.verbose(args,"   Writing "+f_coarse)
        
        nc_coarse = Dataset( f_coarse,'w' )
        nc_coarse.createDimension('time', None )
        nc_coarse.createDimension('y1', coarse.ny)
        nc_coarse.createDimension('x1', coarse.nx)

        coarse.time = nc_coarse.createVariable('time', 'f4', ('time',))
        coarse.y    = nc_coarse.createVariable('y1',   'f4', ('y1',)  )
        coarse.x    = nc_coarse.createVariable('x1',   'f4', ('x1',)  )

        copy_atts(base.y, coarse.y)
        copy_atts(base.x, coarse.x)

        coarse.y[:] = base.y[::skip]
        coarse.x[:] = base.x[::skip]
        coarse.time[0] = 0.

        for var_name, var_data in nc_base.variables.iteritems() : 
            if (var_name != 'x1' and var_name != 'y1' and var_name != 'time'):
                var_coarse = nc_coarse.createVariable(var_name, 'f4', ('time','y1','x1',))
                var_coarse[0,:,:] = var_data[0,::skip,::skip]
                copy_atts(var_data, var_coarse)

        nc_coarse.close()
        # set file permissions
        os.chmod(f_coarse, 0o644)   # uses an Octal number!

        # write config files
        speak.verbose(args,"   Writing the "+str(skip)+" km config file.")
        
        config.write(f_coarse, f_template, coarse, skip)
