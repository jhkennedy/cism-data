import pyproj
import numpy as np
from scipy import interpolate

from util.ncfunc import copy_atts
from util import speak

def get_mcb(args, nc_massCon, nc_base, trans, proj_eigen_gl04c, proj_epsg3413):
    """Get the mass conserving bed data.
    """
    massCon_y = nc_massCon.variables['y']
    massCon_ny = massCon_y[:].shape[0]

    massCon_x = nc_massCon.variables['x']
    massCon_nx = massCon_x[:].shape[0]

    massCon_data = np.ndarray( (massCon_ny,massCon_nx) )
    base_data    = np.ndarray( (trans.ny,trans.nx) )
    temp_data    = np.ndarray( (trans.ny,trans.nx) )

    var_list       = [ 'thickness', 'bed',  'errbed' ]
    rename_massCon = [ 'thk',       'topg', 'topgerr']
    for vv in range(0, len(var_list) ) :
        var = var_list[vv]
        massCon_data[:,:] = 0.
        base_data[:,:]    = 0.
        temp_data[:,:]    = 0.

        massCon_var = nc_massCon.variables[var]
        massCon_data[:,:] = massCon_var[::-1,:] # y fliped when compaired to Bamber

        speak.verbose(args,"   Interpolating "+var_list[vv]+".")
        massCon_to_base = interpolate.RectBivariateSpline( massCon_y[::-1], massCon_x[:], massCon_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster

        for ii in range(0, trans.nx):
            base_data[:,ii] = massCon_to_base.ev(trans.y_grid[:,ii], trans.x_grid[:,ii] )
        
        # This is only needed for data that is actually referenced from the EIGEN-GL04C Geoid -- all topographical 'z' data. Not needed if not 'z' data.
        if (var_list[vv] != 'errbed') :
            temp_x_grid, temp_y_grid, temp_data = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, trans.x_grid.flatten(), trans.y_grid.flatten(), base_data.flatten())
            temp_data = temp_data.reshape((trans.ny,trans.nx))
            base_data[:,:] = temp_data[:,:]
        
        if (var_list[vv] == 'thickness') :
            bad_val = 0
        else :
            bad_val = -9999

        data_min = np.amin( massCon_data[ massCon_data[:,:] != bad_val] )
        data_max = np.amax( massCon_data[ massCon_data[:,:] != bad_val] )

        data_mask = np.ma.masked_outside(base_data, data_min, data_max)
        base_data[ data_mask.mask ] = bad_val

        speak.verbose(args,"   Writing "+rename_massCon[vv]+" to base.")
        base_var = nc_base.createVariable( rename_massCon[vv], 'f4', ('y','x',) )
        base_var[:,:] = base_data[:,:]
        #copy_atts(massCon_var, base_var)
        #NOTE: this data is formated short, while the rest of the data is a float. Does not like _FillValue from short data. 
        atts = massCon_var.ncattrs()
        for ii in range(len(atts)):
            if (atts[ii] != '_FillValue'):
                base_var.setncattr(atts[ii], massCon_var.getncattr(atts[ii]))
            else:
                base_var.setncattr('missing_value', base_data[-1,-1]) # from a known bad point

    # drop temp. variables
    temp_y_grid = None
    temp_x_grid = None
    temp_data   = None
    nc_massCon.close()

