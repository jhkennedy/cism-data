import os
import json
import numpy as np

from netCDF4 import Dataset

from util import speak
from util.ncfunc import copy_atts, get_nc_file
from util.projections import DataGrid
from templates import config

def build_base(f_base, f_epsg, d_meters):
    with open(f_epsg, 'r') as f:
        epsg = json.load(f)

    EPSG = DataGrid()
    EPSG.x = np.arange(epsg['ll'][0], epsg['ur'][0], d_meters)
    EPSG.nx = len(EPSG.x)

    EPSG.y = np.arange(epsg['ll'][1], epsg['ur'][1], d_meters)
    EPSG.ny = len(EPSG.y)

    nc_base = Dataset(f_base,'w', format='NETCDF4')
    nc_base.createDimension('y',EPSG.ny)
    nc_base.createDimension('x',EPSG.nx)
    
    base = DataGrid()
    base.ny = EPSG.ny
    base.nx = EPSG.nx
    
    base.y = nc_base.createVariable('y', 'f4', 'y')
    base.y[:] = EPSG.y[:]
    #FIXME: attributes

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = EPSG.x[:]
    #FIXME: attributes

    # create some grids for interpolation
    base.make_grid()

    return (nc_base, base)


def add_time(args, f_base, f_1km, f_template, f_epsg_shrunk):
    """Add the time dimension to a EPSG:3413 1km dataset and write config files. 

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
    
    base = DataGrid()
    base.y = nc_base.variables['y']
    base.x = nc_base.variables['x']
    
    with open(f_epsg_shrunk, 'r') as f:
        epsg = json.load(f)
    
    idx_left = (np.abs(base.x[:] - epsg['ll'][0])).argmin()
    idx_right = (np.abs(base.x[:] - epsg['ur'][0])).argmin()
    
    idx_lower = (np.abs(base.y[:] - epsg['ll'][1])).argmin()
    idx_upper = (np.abs(base.y[:] - epsg['ur'][1])).argmin()

    y_shrink = [idx_lower,idx_upper+1] #NOTE: python stop exclusive, nco stop inclusive!
    x_shrink = [idx_left, idx_right+1] #NOTE: python stop exclusive, nco stop inclusive!

    if args.extended:
        base.ny = base.y.shape[0]
        base.nx = base.x.shape[0]
    else:
        base.ny = base.y[ y_shrink[0]:y_shrink[1] ].shape[0]
        base.nx = base.x[ x_shrink[0]:x_shrink[1] ].shape[0]

    speak.verbose(args,"   Writing "+f_1km)
    nc_1km = Dataset(f_1km, 'w', format='NETCDF4')
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
            copy_atts(var_data, var_1km)
            vd = var_data[:,:]
            msk = vd < -6000.
            vd[msk] = -9999.
            if args.extended:
                var_1km[0,:,:] = vd[:,:]
            else:
                var_1km[0,:,:] = vd[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]
    #copy_atts(var_data, var_1km)

    nc_base.close()
    nc_1km.close()
    os.chmod(f_1km, 0o644)   # uses an octal number!

    speak.verbose(args,"   Writing the 1km config file.")

    config.write(f_1km, f_template, base, 1)


