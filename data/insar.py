import numpy as np
from scipy import interpolate

from util.ncfunc import copy_atts
from util import speak

def get_velocity(args, nc_insar, nc_base, trans):
    """Get the velocities from the insar data.
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

    nc_insar.close()

