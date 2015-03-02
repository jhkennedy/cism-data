from netCDF4 import Dataset
import numpy as np
"""
A set of useful functions when dealing with netcdf4 data.
"""

def copy_atts(fin, fout) :
    """Copy netCDF attributes.

    This function copies the attributes from one netCDF element to another.
    
    Parameters
    ----------
    fin : 
        Source netCDF element
    fout :
        Target netCDF element

    Examples
    --------
    Copy the attributes from one variable to another.

    >>> old_var = nc_old.variables['old']
    >>> new_var = nc_new.createVariable('new', 'f4', ('y','x',) )
    >>> new_var[:,:] = old_var[:,:]
    >>> copy_atts( old_var,new_var )
    """
    
    # get a list of global attribute names from the incoming file
    atts = fin.ncattrs()
    
    # place those attributes in the outgoing file
    for ii in range(len(atts)) :
        fout.setncattr(atts[ii], fin.getncattr(atts[ii]))


