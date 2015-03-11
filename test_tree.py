import os
import scipy
import pyproj
import numpy as np
from netCDF4 import Dataset
from scipy.spatial import cKDTree

import util.interpolate as interp
from util.ncfunc import copy_atts, get_nc_file
from util.projections import DataGrid


print('Loading Data')

lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'

path_bamber = os.path.dirname(lc_bamber)

nc_bamber = get_nc_file(lc_bamber,'r')
nc_massCon = get_nc_file(lc_massCon,'r')

bamber  = DataGrid()
massCon = DataGrid()

bamber.y = nc_bamber.variables['projection_y_coordinate']
bamber.x = nc_bamber.variables['projection_x_coordinate']
bamber.ny = bamber.y[:].shape[0]
bamber.nx = bamber.x[:].shape[0]
bamber.make_grid()

bamber.yx = np.ndarray( (len( bamber.y_grid.ravel() ),2) )
bamber.yx[:,0] = bamber.y_grid.ravel()
bamber.yx[:,1] = bamber.x_grid.ravel()

massCon.y = nc_massCon.variables['y']
massCon.x = nc_massCon.variables['x']
massCon.ny = massCon.y[:].shape[0]
massCon.nx = massCon.x[:].shape[0]
massCon.make_grid_flip_y()

massCon.yx = np.ndarray( (len( massCon.y_grid.ravel() ),2) )
massCon.yx[:,0] = massCon.y_grid.ravel()
massCon.yx[:,1] = massCon.x_grid.ravel()

print('Creating test data file')
# create some base variables
class test():
    pass

nc_test = Dataset('test.nc','w')
nc_test.createDimension('y',bamber.ny)
nc_test.createDimension('x',bamber.nx)
test.y = nc_test.createVariable('y', 'f4', 'y')
test.y[:] = bamber.y[:]
copy_atts(bamber.y, test.y) #FIXME: units say km, but it's actuall in m

test.x = nc_test.createVariable('x', 'f4', 'x')
test.x[:] = bamber.x[:]
copy_atts(bamber.x, test.x) #FIXME: units say km, but it's actuall in m

print('Building trees')
# ckd trees of data grids
bamber.tree = cKDTree(bamber.yx)
massCon.tree = cKDTree(massCon.yx)

print('Transformations')
# get transformation grids
proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')
proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m')

b_trans = DataGrid()
b_trans.ny = bamber.ny
b_trans.nx = bamber.nx
b_trans.dims = (bamber.ny, bamber.nx)

b_trans.x_grid, b_trans.y_grid = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, bamber.x_grid.flatten(), bamber.y_grid.flatten())
b_trans.y_grid = b_trans.y_grid.reshape((b_trans.ny,b_trans.nx))
b_trans.x_grid = b_trans.x_grid.reshape((b_trans.ny,b_trans.nx))

b_trans.yx = np.ndarray( (len( b_trans.y_grid.ravel() ),2) )
b_trans.yx[:,0] = b_trans.y_grid.ravel()
b_trans.yx[:,1] = b_trans.x_grid.ravel()

b_trans.qd, b_trans.qi = massCon.tree.query(b_trans.yx, k=1)   # nearest neighbor in massCon for transformed bamber grid

print('Bulding data arrays')
# build data arrays
pri_data = np.ma.masked_equal( nc_massCon.variables['bed'][::-1,:], -9999)     # MASSCON GRID
sec_data = np.ma.masked_values( nc_bamber.variables['BedrockElevation'][:,:], -9999.) # BAMBER GRID
new_data = np.ma.array( np.zeros( (bamber.ny,bamber.nx) ), mask=np.zeros( (bamber.ny,bamber.nx) ) )
reason = np.zeros( bamber.dims )

#import matplotlib.pyplot as plt
#plt.imshow(pri_data.filled(-1))
#plt.show()

print('Interpolating...')
# loop through all data points
for ii in range(0, bamber.ny):
    for jj in range(0, bamber.nx):
        # interpolate bamber point on massCon grid
        indx = np.ravel_multi_index( (ii,jj),b_trans.dims )

        nn_mc = b_trans.qi[indx]
        nn_ii, nn_jj = np.unravel_index(nn_mc, massCon.dims)
        
        y_s = -1
        x_s = -1
        
        # make sure inside priority grid (massCon)
        if (b_trans.y_grid[ii,jj] < massCon.y_grid[0,0] or b_trans.y_grid[ii,jj] > massCon.y_grid[-1,0]):
            # outside y range
            if sec_data.mask[ii,jj]:
                new_data.mask[ii,jj] = True
                reason[ii,jj] = 1
                continue
            else:
                tx, ty, td = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, bamber.x_grid[ii,jj], bamber.y_grid[ii,jj], sec_data[ii,jj])
                new_data[ii,jj] = td
                reason[ii,jj] = -1
                continue
        
        if (b_trans.x_grid[ii,jj] < massCon.x_grid[0,0] or b_trans.x_grid[ii,jj] > massCon.x_grid[0,-1]):
            # outside x range 
            if sec_data.mask[ii,jj]:
                new_data.mask[ii,jj] = True
                reason[ii,jj] = 2
                continue
            else:
                tx, ty, td = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, bamber.x_grid[ii,jj], bamber.y_grid[ii,jj], sec_data[ii,jj])
                new_data[ii,jj] = td
                reason[ii,jj] = -1
                continue

        # find quadrent
        if b_trans.y_grid[ii,jj] >= massCon.y_grid[nn_ii,nn_jj]:
            y_s = +1
        if b_trans.x_grid[ii,jj] >= massCon.x_grid[nn_ii,nn_jj]:
            x_s = +1

        # check for missing priority data!
        missing_points, interp_dict = interp.check_missing(pri_data, (nn_ii, nn_jj), y_s, x_s)

        # get secondary data!
        if not missing_points :
            pass

        elif not sec_data.mask[ii,jj] :
            tx, ty, td = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, bamber.x_grid[ii,jj], bamber.y_grid[ii,jj], sec_data[ii,jj])
            if len(missing_points) <= 3 :
                for point in missing_points:
                    # use secondary data at (ii,jj) for missing points, but keep same interp weight!
                    interp_dict[point] = td
                reason[ii,jj] = -2
            
            else:
                new_data[ii,jj] = td
                reason[ii,jj] = -1

        else:
            new_data.mask[ii,jj] = True
            reason[ii,jj] = 3
            continue

        # interpolate!
        alpha = ( b_trans.y_grid[ii,jj] - massCon.y_grid[nn_ii,nn_jj] )/(massCon.dy*y_s)
        beta  = ( b_trans.x_grid[ii,jj] - massCon.x_grid[nn_ii,nn_jj] )/(massCon.dx*x_s)

        w = interp.linear_weights(alpha, beta)

        new_data[ii,jj] = interp_dict[ (nn_ii,    nn_jj    ) ]*w[0] \
                         +interp_dict[ (nn_ii,    nn_jj+x_s) ]*w[1] \
                         +interp_dict[ (nn_ii+y_s,nn_jj+x_s) ]*w[2] \
                         +interp_dict[ (nn_ii+y_s,nn_jj    ) ]*w[3]

import matplotlib.pyplot as plt
plt.imshow(reason)
plt.show()

print('Writing data')
test_var = nc_test.createVariable( 'interpData', 'f4', ('y','x',) )
test_var[:,:] = new_data[:,:]

test_reason = nc_test.createVariable('reason','i',('y','x',))
test_reason[:,:] = reason[:,:]
nc_test.close()
print('And Done!')
