import os

from netCDF4 import Dataset

from util import speak
from util.ncfunc import copy_atts, get_nc_file
from util.projections import DataGrid
from templates import config

def build_base (nc_bamber, nc_base, base):
    """Build the Bamber base grid.
    """
    bamber_y = nc_bamber.variables['projection_y_coordinate']
    bamber_ny = bamber_y[:].shape[0] # number of y points for 1km grid

    bamber_x = nc_bamber.variables['projection_x_coordinate']
    bamber_nx = bamber_x[:].shape[0] # number of x points for 1km grid

    # make bamber 1km grid for base
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

    return base

def add_time (args, f_base, f_1km, f_template):
    """Add the time dimention to a Bamber 1km DEM dataset and write config files. 
    """
    # shrink dataset to the ice sheet
    nc_base = Dataset(f_base,'r')
    
    y_shrink = [100,2900+1] #NOTE: python stop exclusive, nco stop inclusive!
    x_shrink = [500,2000+1] #NOTE: python stop exclusive, nco stop inclusive!

    base = DataGrid()
    base.y = nc_base.variables['y']
    base.ny = base.y[ y_shrink[0]:y_shrink[1] ].shape[0]
    base.x = nc_base.variables['x']
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

    y1[:] = base.y[ y_shrink[0]:y_shrink[1] ]
    x1[:] = base.x[ x_shrink[0]:x_shrink[1] ]
    time[0] = 0.

    for var_name, var_data in nc_base.variables.iteritems() : 
        if (var_name != 'x' and var_name != 'y'):
            var_1km = nc_1km.createVariable(var_name, 'f4', ('time','y1','x1',))
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
