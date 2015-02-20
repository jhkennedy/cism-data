# a set of useful functions when dealing with netcdf4 Data

from netCDF4 import Dataset
import numpy as np

# copy nc global attributes
def copy_atts(fin, fout) :
    # get a list of global attribute names from the incoming file
    atts = fin.ncattrs()
    # place those attributes in the outgoing file
    for ii in range(len(atts)) :
        fout.setncattr(atts[ii], fin.getncattr(atts[ii]))


