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
    base.y.long_name = 'y-coordinate in projection'
    base.y.standard_name = 'projection_y_coordinate'
    base.y.axis = 'Y'
    base.y.units = 'meters'

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = EPSG.x[:]
    base.x.long_name = 'x-coordinate in projection'
    base.x.standard_name = 'projection_x_coordinate'
    base.x.axis = 'X'
    base.x.units = 'meters'
    
    base.proj = nc_base.createVariable('epsg_3413', 'b')
    base.proj.grid_mapping_name = 'polar_stereographic' 
    base.proj.latitude_of_projection_origin = '+90' 
    base.proj.straight_vertical_longitude_from_pole = '-45' 
    base.proj.standard_parallel = '70' 
    base.proj.proj_scale_factor = '1' 
    base.proj.false_easting = '0' 
    base.proj.false_northing = '0' 
    base.proj.ellipsoid = 'WGS84' 
    base.proj.datum = 'WGS84' 
    base.proj.units = 'meters' 
    base.proj.proj4_string = '+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' 

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

    base.ny = base.y[ y_shrink[0]:y_shrink[1] ].shape[0]
    base.nx = base.x[ x_shrink[0]:x_shrink[1] ].shape[0]

    speak.verbose(args,"   Writing "+f_1km)
    nc_1km = Dataset(f_1km, 'w', format='NETCDF4')
    nc_1km.createDimension('time', None )
    nc_1km.createDimension('y1', base.ny)
    nc_1km.createDimension('x1', base.nx)

    time = nc_1km.createVariable('time', 'f4', ('time',))
    time.long_name = 'time'
    time.units = 'common_years since 2009-01-01 00:00:00'
    time.calendar = '365_day'
    time.comment = "The initial time here is an estimate of the nominal date for Joughin's 2015 "+ \
                   "InSAR data. Because this is a synthesis of datasets across many time periods, "+ \
                   "the inital date is inherently fuzzy and should be changed to suit your purposes."
    
    y1   = nc_1km.createVariable('y1',   'f4', ('y1',)  )
    x1   = nc_1km.createVariable('x1',   'f4', ('x1',)  )
    copy_atts(base.y, y1)
    copy_atts(base.x, x1)

    y1[:] = base.y[ y_shrink[0]:y_shrink[1] ]
    x1[:] = base.x[ x_shrink[0]:x_shrink[1] ]
    time[0] = 0.

    base.proj = nc_base.variables['epsg_3413']
    proj = nc_1km.createVariable('epsg_3413', 'b')
    copy_atts(base.proj, proj)
    
    for var_name, var_data in nc_base.variables.iteritems() : 
        if var_name not in  ['x', 'y', 'lat', 'lon', 'epsg_3413']:
            var_1km = nc_1km.createVariable(var_name, 'f4', ('time','y1','x1',))
            copy_atts(var_data, var_1km)
            var_1km[0,:,:] = var_data[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]
        
        elif var_name not in ['x', 'y', 'epsg_3413']:
            var_1km = nc_1km.createVariable(var_name, 'f4', ('y1','x1',))
            copy_atts(var_data, var_1km)
            var_1km[:,:] = var_data[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]

    nc_base.close()
    nc_1km.close()
    os.chmod(f_1km, 0o644)   # uses an octal number!

    speak.verbose(args,"   Writing the 1km config file.")

    config.write(f_1km, f_template, base, 1)


