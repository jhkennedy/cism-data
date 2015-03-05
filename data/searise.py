import numpy as np

from util.ncfunc import copy_atts
from util import speak

def get_bheatflx_artm(args, nc_seaRise, nc_base, base) :
    """Get bheatflx and artm from the sea rise data.
    """
    seaRise_y = nc_seaRise.variables['y']
    seaRise_ny = seaRise_y[:].shape[0]

    seaRise_x = nc_seaRise.variables['x']
    seaRise_nx = seaRise_x[:].shape[0]

    # to convert to Bamber 1km grid
    seaRise_data = np.ndarray( (base.ny, base.nx) )
    seaRise_y_equal_base = 100
    seaRise_x_equal_base = 500

    # get basal heat flux
    #--------------------
    seaRise_data[:,:] = 0.
    seaRise_bheatflx = nc_seaRise.variables['bheatflx']
    seaRise_data[ seaRise_y_equal_base:seaRise_y_equal_base+seaRise_ny , seaRise_x_equal_base:seaRise_x_equal_base+seaRise_nx ] = -seaRise_bheatflx[0,:,:] # invert sign!

    speak.verbose(args,"   Writing bheatflx to base.")
    base_bheatflx = nc_base.createVariable('bheatflx', 'f4', ('y','x',) )
    base_bheatflx[:,:] = seaRise_data[:,:] 
    copy_atts(seaRise_bheatflx, base_bheatflx)

    # get annual mean air temperature (2m)
    #-------------------------------------
    seaRise_data[:,:] = 0.
    seaRise_presartm = nc_seaRise.variables['presartm']
    seaRise_data[ seaRise_y_equal_base:seaRise_y_equal_base+seaRise_ny , seaRise_x_equal_base:seaRise_x_equal_base+seaRise_nx ] = seaRise_presartm[0,:,:]

    speak.verbose(args,"   Writing artm to base.")
    base_artm = nc_base.createVariable('artm', 'f4', ('y','x',) )
    base_artm[:,:] = seaRise_data[:,:] 
    copy_atts(seaRise_presartm, base_artm)

