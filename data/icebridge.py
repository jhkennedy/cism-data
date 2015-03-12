import pyproj
import scipy
import numpy as np

from util import speak
from util.ncfunc import copy_atts, copy_atts_bad_fill
from util.projections import DataGrid
import util.interpolate as interp

def get_mcb(args, nc_massCon, nc_bamber, nc_base, base, trans, proj_eigen_gl04c, proj_epsg3413):
    """Get the mass conserving bed data.
    """
    massCon = DataGrid()
    
    massCon.y = nc_massCon.variables['y']
    massCon.x = nc_massCon.variables['x']
    massCon.ny = massCon.y[:].shape[0]
    massCon.nx = massCon.x[:].shape[0]
    massCon.make_grid_flip_y()

    massCon.yx = np.ndarray( (len(massCon.y_grid.ravel()),2) )
    massCon.yx[:,0] = massCon.y_grid.ravel()
    massCon.yx[:,1] = massCon.x_grid.ravel()
    
    massCon.tree = scipy.spatial.cKDTree(massCon.yx)

    trans.yx = np.ndarray( (len(trans.y_grid.ravel()),2) )
    trans.yx[:,0] = trans.y_grid.ravel()
    trans.yx[:,1] = trans.x_grid.ravel()

    trans.qd, trans.qi = massCon.tree.query(trans.yx, k=1)   # nearest neighbor in massCon for transformed base grid

    massCon.thickness = nc_massCon.variables['thickness']
    massCon.thk = np.ndarray(massCon.dims)
    massCon.thk[:,:] = massCon.thickness[::-1,:] # y fliped when compaired to Bamber

    speak.verbose(args,"   Interpolating thickness.")
    massCon_to_base = scipy.interpolate.RectBivariateSpline( massCon.y[::-1], massCon.x[:], massCon.thk, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
    trans.thk = np.zeros( trans.dims )
    for ii in range(0, base.nx):
        trans.thk[:,ii] = massCon_to_base.ev(trans.y_grid[:,ii], trans.x_grid[:,ii] )
    
   
    speak.verbose(args,"   Writing thk to base.")
    base.thk = nc_base.createVariable('thk', 'f4', ('y','x',) )
    base.thk[:,:] = trans.thk[:,:]
    copy_atts_bad_fill(massCon.thickness, base.thk, -9999.)

    speak.verbose(args,"   Interpolating, with priority, topg and topgerr.")
    speak.verbose(args,"      Primary Data [IceBridge]:  bed and errbed.")
    speak.verbose(args,"      Secondary Data [BamberDEM]:  BedrockElevation and BedrockError.")

    pri_data = np.ma.masked_equal( nc_massCon.variables['bed'][::-1,:] , -9999)
    sec_data = np.ma.masked_values( nc_bamber.variables['BedrockElevation'][:,:], -9999.)
    new_data = np.ma.array(np.zeros(base.dims), mask=np.zeros(base.dims))

    pri_err = np.ma.masked_equal( nc_massCon.variables['errbed'][::-1,:] , -9999)
    sec_err = np.ma.masked_values( nc_bamber.variables['BedrockError'][:,:], -9999.)
    new_err = np.ma.array(np.zeros(base.dims), mask=np.zeros(base.dims))

    for ii in range(0, base.ny):
        for jj in range(0, base.nx):
            # make sure inside priority grid (massCon)
            if (trans.y_grid[ii,jj] < massCon.y_grid[0,0] or trans.y_grid[ii,jj] > massCon.y_grid[-1,0]):
                # outside y range
                if sec_data.mask[ii,jj]:
                    new_data.mask[ii,jj] = True
                    new_err.mask[ii,jj] = True
                    continue
                else:
                    tx, ty, td = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, base.x_grid[ii,jj], base.y_grid[ii,jj], sec_data[ii,jj])
                    new_data[ii,jj] = td
                    new_err[ii,jj] = sec_err[ii,jj]
                    continue
            
            if (trans.x_grid[ii,jj] < massCon.x_grid[0,0] or trans.x_grid[ii,jj] > massCon.x_grid[0,-1]):
                # outside x range 
                if sec_data.mask[ii,jj]:
                    new_data.mask[ii,jj] = True
                    new_err.mask[ii,jj] = True
                    continue
                else:
                    tx, ty, td = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, base.x_grid[ii,jj], base.y_grid[ii,jj], sec_data[ii,jj])
                    new_data[ii,jj] = td
                    new_err[ii,jj] = sec_err[ii,jj]
                    continue

            indx = np.ravel_multi_index( (ii,jj), base.dims)
            
            # nearest neighbor indices
            nn_ii, nn_jj = np.unravel_index( trans.qi[indx], massCon.dims)

            # to find nearest neighbors quadrent
            i_s = -1
            j_s = -1

            # find quadrent
            if trans.y_grid[ii,jj] >= massCon.y_grid[nn_ii,nn_jj]:
                i_s = +1
            if trans.x_grid[ii,jj] >= massCon.x_grid[nn_ii,nn_jj]:
                j_s = +1

            # check for missing priority data!
            missing_points, interp_dict = interp.check_missing(pri_data, (nn_ii, nn_jj), i_s, j_s)
            missing_err_pts, err_dict = interp.check_missing(pri_err, (nn_ii, nn_jj), i_s, j_s)

            # get secondary data!
            if not missing_points :
                pass

            elif not sec_data.mask[ii,jj] :
                tx, ty, td = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, base.x_grid[ii,jj], base.y_grid[ii,jj], sec_data[ii,jj])
                if len(missing_points) <= 3 :
                    for point in missing_points:
                        # use secondary data at (ii,jj) for missing points, but keep same interp weight!
                        interp_dict[point] = td
                        err_dict[point] = sec_err[ii,jj]
                
                else:
                    new_data[ii,jj] = td
                    new_err[ii,jj] = sec_err[ii,jj]
                    continue

            else:
                new_data.mask[ii,jj] = True
                new_err.mask[ii,jj] = True
                continue

            # interpolate!
            alpha = ( trans.y_grid[ii,jj] - massCon.y_grid[nn_ii,nn_jj] )/(massCon.dy*i_s)
            beta  = ( trans.x_grid[ii,jj] - massCon.x_grid[nn_ii,nn_jj] )/(massCon.dx*j_s)

            w = interp.linear_weights(alpha, beta)

            new_data[ii,jj] = interp_dict[ (nn_ii,    nn_jj    ) ]*w[0] \
                             +interp_dict[ (nn_ii,    nn_jj+j_s) ]*w[1] \
                             +interp_dict[ (nn_ii+i_s,nn_jj+j_s) ]*w[2] \
                             +interp_dict[ (nn_ii+i_s,nn_jj    ) ]*w[3]


            new_err[ii,jj] = err_dict[ (nn_ii,    nn_jj    ) ]*w[0] \
                            +err_dict[ (nn_ii,    nn_jj+j_s) ]*w[1] \
                            +err_dict[ (nn_ii+i_s,nn_jj+j_s) ]*w[2] \
                            +err_dict[ (nn_ii+i_s,nn_jj    ) ]*w[3]

            missing_points = None
            interp_dict = None
            err_dict = None
    # now transform new data back to bamber grid. 
    temp_x_grid, temp_y_grid, temp_data = pyproj.transform(proj_epsg3413, proj_eigen_gl04c, trans.x_grid.flatten(), trans.y_grid.flatten(), new_data.flatten())
    new_data[:,:] = temp_data.reshape( base.dims )[:,:]

    speak.verbose(args,"   Writing topg topgerr to base.")
    base.topg = nc_base.createVariable('topg', 'f4', ('y','x',) )
    base.topg[:,:] = new_data.filled(-9999.)[:,:]
    copy_atts_bad_fill(nc_massCon.variables['bed'], base.topg, -9999.)

    base.topgerr = nc_base.createVariable('topgerr', 'f4', ('y','x',) )
    base.topgerr[:,:] = new_err.filled(-9999.)[:,:]
    copy_atts_bad_fill(nc_massCon.variables['errbed'], base.topgerr, -9999.)
    

