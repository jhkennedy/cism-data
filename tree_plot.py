import numpy as np
from util.ncfunc import copy_atts, get_nc_file
from util.projections import DataGrid


nc_test = get_nc_file('test.nc','r')
test = DataGrid()

test.reason = nc_test.variables['reason']
test.dims = test.reason[:,:].shape

reason = np.ndarray(test.dims)
reason[:,:] = test.reason[:,:]

test.img = np.ones( (test.dims[0], test.dims[1], 3) )
test.img *= 255

# use secondary as interp data! 
test.img[reason == -3, : ] = [255,0,0]
# use secondary as new data! 
test.img[reason == -1, : ] = [102,255,102]
test.img[reason == -2, : ] = [102,255,102]
# outside y range and no secondary data
test.img[reason == 1, : ] = [255,102,255]
# outside x range and no secondary data 
test.img[reason == 2, : ] = [178,102,255]
# inside but no secondary data
test.img[reason == 3, : ] = [255,102,178]



test.img /= 255

import matplotlib.pyplot as plt
plt.imshow(test.img[::-1,:,:])
plt.show()

nc_test.close
