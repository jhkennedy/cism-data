import numpy as np

from util.ncfunc import copy_atts
from util import speak

def get_acab(args, nc_racmo2p0, nc_base, base):
    """Get acab from the RACMO2.0 data.
    """
    racmo2p0_data = np.ndarray( (base.ny, base.nx) )

    racmo2p0_data[:,:] = 0.
    racmo2p0_smb = nc_racmo2p0.variables['smb']
    racmo2p0_data[:,:] = racmo2p0_smb[:,::-1].transpose() / 910.
    racmo2p0_data = np.ma.masked_invalid(racmo2p0_data) # find invalid data and create a mask
    racmo2p0_data = racmo2p0_data.filled(0.)            # fill invalid data with zeros

    speak.verbose(args,"   Writing acab to base.")
    base_acab = nc_base.createVariable( 'acab', 'f4', ('y','x',) )
    base_acab[:,:] = racmo2p0_data[:,:]
    copy_atts(racmo2p0_smb, base_acab) #FIXME: check atribute units -- divided by 910 earlier

