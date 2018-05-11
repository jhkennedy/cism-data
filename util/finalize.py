import os
import json
import numpy as np

from netCDF4 import Dataset

from util import speak
from util.ncfunc import copy_atts
from util.projections import DataGrid
from templates import config

def add_time_and_shrink(args, proj_name, f_base, f_1km, f_template, f_shrink):
    """Add the time dimension to a EPSG:3413 1km dataset and write config files. 

    This function opens the base dataset, creates a new 1km dataset with the time
    dimension that has been shrunken to just around the ice sheet, and creates a 
    CISM config file. NetCDF attributes are preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    proj_name :
        Name of the projection the dataset was created with.
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
    
    with open(f_shrink, 'r') as f:
        shrink = json.load(f)
    
    idx_left = (np.abs(base.x[:] - shrink['ll'][0])).argmin()
    idx_right = (np.abs(base.x[:] - shrink['ur'][0])).argmin()
    
    idx_lower = (np.abs(base.y[:] - shrink['ll'][1])).argmin()
    idx_upper = (np.abs(base.y[:] - shrink['ur'][1])).argmin()

    y_shrink = [idx_lower,idx_upper+1]  # NOTE: python stop exclusive, nco stop inclusive!
    x_shrink = [idx_left, idx_right+1]  # NOTE: python stop exclusive, nco stop inclusive!

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

    if nc_base.variables.has_key(proj_name):
        base.proj = nc_base.variables[proj_name]
        proj = nc_1km.createVariable(proj_name, 'b')
        copy_atts(base.proj, proj)
    
    for var_name, var_data in nc_base.variables.iteritems() : 
        if var_name not in  ['x', 'y', 'lat', 'lon', proj_name]:
            var_1km = nc_1km.createVariable(var_name, 'f4', ('time','y1','x1',))
            copy_atts(var_data, var_1km)
            var_1km[0,:,:] = var_data[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]
        
        elif var_name not in ['x', 'y', 'epsg_3413', proj_name]:
            var_1km = nc_1km.createVariable(var_name, 'f4', ('y1','x1',))
            copy_atts(var_data, var_1km)
            var_1km[:,:] = var_data[ y_shrink[0]:y_shrink[1] , x_shrink[0]:x_shrink[1] ]

    nc_base.close()
    nc_1km.close()
    os.chmod(f_1km, 0o644)   # uses an octal number!

    speak.verbose(args,"   Writing the 1km config file.")

    config.write(f_1km, f_template, base, 1)


def coarsen(args, proj_name, f_base, f_template, coarse_list):
    """Coarsen an base dataset.

    This function opens the 1 km dataset and creates coarser resolution datasets.
    All NetCDF data attributes are preserved.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    proj_name :
        Name of the projection the dataset was created with.
    f_base :
        Filename for the base dataset to coarsen.
    f_template :
        Filename for the template used to create the CISM config files.
    coarse_list :
        List of resolutions to coarsen to. 
    """
    nc_base = Dataset(f_base,'r' )
    
    base = DataGrid()
    base.time = nc_base.variables['time']
    base.y = nc_base.variables['y1']
    base.ny = base.y[:].shape[0]
    base.x = nc_base.variables['x1']
    base.nx = base.x[:].shape[0]

    base.proj = nc_base.variables[proj_name]
    
    for skip in coarse_list :
        coarse = DataGrid()

        coarse.ny = base.y[::skip].shape[0]
        coarse.nx = base.x[::skip].shape[0]
       
        idx = f_base.find('km')
        f_coarse = f_base[:idx-1]+str(skip)+f_base[idx:]

        speak.verbose(args,"   Writing "+f_coarse)
        
        nc_coarse = Dataset( f_coarse,'w', format='NETCDF4' )
        copy_atts(nc_base, nc_coarse)

        nc_coarse.createDimension('time', None )
        nc_coarse.createDimension('y1', coarse.ny)
        nc_coarse.createDimension('x1', coarse.nx)

        coarse.time = nc_coarse.createVariable('time', 'f4', ('time',))
        coarse.y    = nc_coarse.createVariable('y1',   'f4', ('y1',)  )
        coarse.x    = nc_coarse.createVariable('x1',   'f4', ('x1',)  )

        copy_atts(base.time, coarse.time)
        copy_atts(base.y, coarse.y)
        copy_atts(base.x, coarse.x)

        coarse.y[:] = base.y[::skip]
        coarse.x[:] = base.x[::skip]
        coarse.time[0] = 0.

        coarse.proj = nc_coarse.createVariable(proj_name, 'c')
        copy_atts(base.proj, coarse.proj)

        for var_name, var_data in nc_base.variables.iteritems():
            if var_name not in  ['time', 'x1', 'y1', 'lat', 'lon', proj_name]:
                var_coarse = nc_coarse.createVariable(var_name, 'f4', ('time','y1','x1',))
                var_coarse[0,:,:] = var_data[0,::skip,::skip]
                copy_atts(var_data, var_coarse)
            
            elif var_name not in ['time', 'x1', 'y1', proj_name]:
                var_coarse = nc_coarse.createVariable(var_name, 'f4', ('y1','x1',))
                var_coarse[:,:] = var_data[::skip,::skip]
                copy_atts(var_data, var_coarse)

        nc_coarse.close()
        # set file permissions
        os.chmod(f_coarse, 0o644)   # uses an Octal number!

        # write config files
        speak.verbose(args,"   Writing the "+str(skip)+" km config file.")
        
        # config.write(f_coarse, f_template, coarse, skip)
