"""
data.icebridge : IceBridge BedMachine Greenland data import module.

This module provides functions to import data from the IceBrige BedMachine 
Greenland dataset. 

Functions list:
    * get_mcb(args, nc_massCon, nc_bamber, nc_base, base, trans, 
              proj_eigen_gl04c, proj_epsg3413)

Notes
-----
This dataset contains bed topography beneath the Greenland ice sheet based on 
mass conservation derived from airborne radar tracks and satellite radar. The 
data set also includes surface and ice thickness measurements.

More information can be found at:

http://nsidc.org/data/docs/daac/icebridge/idbmg4/index.html

The data uses the ESPG:3413 projection:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = 70 degrees
    * Latitude of projection origin = 90 degrees
    * Central meridian = -45 degrees
    * false eastings = 0
    * flase northings = 0
    * 150 m postings with
        + lower-left corner y,x: -3349425,-637925 (m) 
        + upper-right corner y,x: -657675, 864625 (m)

The dataset is on a 10018 x 17946 grid and the data is formated as a short number.

References
----------
M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour, Deeply incised 
submarine glacial valleys beneath the Greenland Ice Sheet, Nat. Geosci., 7, 
418-422, 2014, doi:10.1038/ngeo2167, 
http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html
"""
import sys
import scipy
import pyproj
import numpy as np
from shapely.geometry import Point, shape 

from util import speak
from util.ncfunc import copy_atts, copy_atts_bad_fill, copy_atts_add_fill
from util import projections 
import util.interpolate as interp

def mcb_epsg3413(args, nc_massCon, nc_bamber, nc_base, base, proj_epsg3413, proj_eigen_gl04c):
    """The mass conserving bed data on CISM's ESPG:3413 grid.

    This function pulls in the `thickness` variable from the mass conserving
    bed dataset, interpolates it to CISM's EPSG:3413 grid, and writes it to the
    base dataset as `thk`. NetCDF attributes are mostly preserved, but the data
    is changed from type short to type float. 
    """
    massCon = projections.DataGrid()
    massCon.y = nc_massCon.variables['y']
    massCon.x = nc_massCon.variables['x']
    massCon.ny = massCon.y[:].shape[0]
    massCon.nx = massCon.x[:].shape[0]
    massCon.make_grid_flip_y()

    massCon.thickness = nc_massCon.variables['thickness']
    massCon.thk = np.ndarray(massCon.dims)
    massCon.thk[:,:] = massCon.thickness[::-1,:] # y fliped

    bamber = projections.DataGrid()
    bamber.y = nc_bamber.variables['projection_y_coordinate']
    bamber.x = nc_bamber.variables['projection_x_coordinate']
    bamber.ny = bamber.y[:].shape[0]
    bamber.nx = bamber.x[:].shape[0] 
    bamber.make_grid()

    speak.verbose(args,"   Interpolating thickness and writing to base.")
    sys.stdout.write(  "   [%-60s] %d%%" % ('='*0, 0.))
    sys.stdout.flush()
    massCon_to_base = scipy.interpolate.RectBivariateSpline( massCon.y[::-1], massCon.x[:], massCon.thk, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
    base.thk = nc_base.createVariable('thk', 'f4', ('y','x',) )
    base.thk[:] = np.zeros( base.dims )
    for ii in range(0, base.nx):
        ctr = (ii*60)/base.nx
        sys.stdout.write("\r   [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
        sys.stdout.flush()
        base.thk[:,ii] = massCon_to_base.ev(base.y_grid[:,ii], base.x_grid[:,ii] )
    sys.stdout.write("\r   [%-60s] %d%%\n" % ('='*60, 100.))
    sys.stdout.flush()
    copy_atts_bad_fill(massCon.thickness, base.thk, -9999.)
    base.thk.grid_mapping = 'epsg_3413'
    base.thk.coordinates = 'lon lat'
    base.thk.reference = 'M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour, Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet, Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167, http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html'


    speak.verbose(args,"   Interpolating, with priority, topg and topgerr.")
    speak.verbose(args,"      Primary Data [IceBridge]:  bed and errbed.")
    speak.verbose(args,"      Secondary Data [BamberDEM]:  BedrockElevation and BedrockError.")

    #base_vars = {'topg':['bed','BedrockElevation']} #NOTE: topgerr accounts for ~10 min. of the runtime.
    base_vars = {'topg':['bed','BedrockElevation'],
                 'topgerr':['errbed', 'BedrockError']}
    for var, var_list in base_vars.iteritems():
        speak.verbose(args,  '\n      Begin '+var+':')
        pri_data = np.ma.masked_equal( nc_massCon.variables[var_list[0]][::-1,:] , -9999)
        sec_data = np.ma.masked_values( nc_bamber.variables[var_list[1]][:,:], -9999.)
        
        pri_range = [pri_data.min(), pri_data.max()]
        rng = [sec_data.min(), sec_data.max()]
        if pri_range[0] < rng[0]:
            rng[0] = pri_range[0]
        if pri_range[1] > rng[1]:
            rng[1] = pri_range[1]

        # fill in missing data values in the IceBridge data with the Bamber DEM
        speak.verbose(args,"         Combining IceBridge and Bamber DEMs.")
        sys.stdout.write(  "         [%-60s] %d%%" % ('='*0, 0.))
        sys.stdout.flush()
        massCon2bamber = projections.transform(massCon, proj_epsg3413, proj_eigen_gl04c)
        sec_interp = scipy.interpolate.RectBivariateSpline( bamber.y[::], bamber.x[:], sec_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
        massCon_bamber = np.zeros( massCon.dims )
        for ii in range(0, massCon.nx):
            ctr = (ii*60)/massCon.nx
            sys.stdout.write("\r         [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
            sys.stdout.flush()
            massCon_bamber[:,ii] = sec_interp.ev(massCon2bamber.y_grid[:,ii], massCon2bamber.x_grid[:,ii] )
        sys.stdout.write("\r         [%-60s] %d%%\n" % ('='*60, 100.))
        sys.stdout.flush()
        pri_data.unshare_mask()
        pri_data[pri_data.mask] = massCon_bamber[pri_data.mask]

        # interpolate the Bamber DEM data to the base grid
        speak.verbose(args,"         Interpolating extended Bamber DEM to base grid.")
        sys.stdout.write(  "         [%-60s] %d%%" % ('='*0, 0.))
        sys.stdout.flush()
        base2bamber = projections.transform(base, proj_epsg3413, proj_eigen_gl04c)
        base_bamber = np.zeros( base.dims )
        for ii in range(0, base.nx):
            ctr = (ii*60)/base.nx
            sys.stdout.write("\r         [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
            sys.stdout.flush()
            base_bamber[:,ii] = sec_interp.ev(base2bamber.y_grid[:,ii], base2bamber.x_grid[:,ii] )
        sys.stdout.write("\r         [%-60s] %d%%\n" % ('='*60, 100.))
        sys.stdout.flush()

        # interpolate the filled IceBridge data to the base grid
        speak.verbose(args,"         Interpolating combined dataset to base grid.")
        sys.stdout.write(  "         [%-60s] %d%%" % ('='*0, 0.))
        sys.stdout.flush()
        pri_interp = scipy.interpolate.RectBivariateSpline( massCon.y[::-1], massCon.x[:], pri_data, kx=1, ky=1, s=0) # regular 2d linear interp. but faster
        base_mcb = np.zeros( base.dims )
        for ii in range(0, base.nx):
            ctr = (ii*60)/base.nx
            sys.stdout.write("\r         [%-60s] %d%%" % ('='*ctr, ctr/60.*100.))
            sys.stdout.flush()
            base_mcb[:,ii] = pri_interp.ev(base.y_grid[:,ii], base.x_grid[:,ii] )
        sys.stdout.write("\r         [%-60s] %d%%\n" % ('='*60, 100.))
        sys.stdout.flush()

        #NOTE: RectBivariateSpline extrapolates data outside the convex hull
        #      (constant value), so, we need to create a mask for values outisde
        #      the IceBridge convex hull...
        base_y_mask = np.ma.masked_outside(base.y_grid,massCon.y[0], massCon.y[-1])
        base_x_mask = np.ma.masked_outside(base.x_grid,massCon.x[0], massCon.x[-1])
        base_mcb_masked = np.ma.masked_array(base_mcb, mask=np.logical_or(base_y_mask.mask,base_x_mask.mask))

        base_bamber[~base_mcb_masked.mask] = base_mcb_masked[~base_mcb_masked.mask]

        bamx = [bamber.x[0], bamber.x[-1], bamber.x[-1], bamber.x[0], bamber.x[0]]
        bamy = [bamber.y[0], bamber.y[0], bamber.y[-1], bamber.y[-1], bamber.y[0]]
        bam_shape = shape({'type':'polygon', 'coordinates':[zip(bamx,bamy)]})
        base_pts_in_bamber = zip(base2bamber.x_grid.flatten(), base2bamber.y_grid.flatten())
        msk = np.ones(base2bamber.x_grid.size, dtype=bool)
        for ii in range(msk.size):
            msk[ii] = Point(base_pts_in_bamber[ii]).within(bam_shape)
        msk.shape = base2bamber.x_grid.shape

        #from pprint import pprint as pp
        #print(msk.size)
        #pp(msk)

        base_bamber[~msk] = -9999.


        #NOTE: Make sure all values fall within a reasonable range as 
        #      RectBivariateSpine interps using the missing values
        base_bamber[base_bamber < rng[0]] = -9999.
        base_bamber[base_bamber > rng[1]] = -9999.

        base.var = nc_base.createVariable(var, 'f4', ('y','x',) )
        base.var[:] = base_bamber[:]  
        copy_atts_bad_fill(nc_massCon.variables[var_list[0]], base.var, -9999.)
        base.var.grid_mapping = 'epsg_3413'
        base.var.coordinates = 'lon lat'
        base.var.reference = 'M. Morlighem, E. Rignot, J. Mouginot, H. Seroussi and E. Larour, Deeply incised submarine glacial valleys beneath the Greenland Ice Sheet, Nat. Geosci., 7, 418-422, 2014, doi:10.1038/ngeo2167, http://www.nature.com/ngeo/journal/vaop/ncurrent/full/ngeo2167.html'
        


def mcb_bamber(args, nc_massCon, nc_bamber, nc_base, base, trans, proj_eigen_gl04c, proj_epsg3413):
    """The mass conserving bed data in the bamber projection.

    This function pulls in the `thickness` variable from the mass conserving
    bed dataset, interpolates it to the Bamber DEM, and writes it to the base 
    dataset as `thk`. NetCDF attributes are mostly preserved, but the data is
    changed from type short to type float. 

    This function then pulls in the `bed` and `errbed` variables from the mass
    conserving bed dataset, and the `bedrockElevation` and `bedrockError`
    variables from the Bamber dataset. Prioritizing the mass conserving bed data,
    the bedrock elevation and error is interpolated to the Bamber DEM and written
    as `topg` and `topgerr` respectively. Data attributes are taken from the mass
    conserving bed data, and the data is written as a float. 

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments. 
    nc_massCon :
        An opened netCDF Dataset containing the Sea Rise data.
    nc_bamber :
        An opened netCDF Dataset containing the Sea Rise data.
    nc_base :
        The created netCDF Dataset that will contain the base data.
    base :
        A DataGrid() class instance that holds the base data grid information.
    trans :
        A DataGrid() class instance that holds the base data grid transformed 
        to EPSG:3413.
    proj_eigen_gl04c :
        A proj class instance that holds the Bamber DEM projection.
    proj_epsg3413 :
        A proj class instance that holds the EPSG:3413 projection.
    """
    massCon = projections.DataGrid()
    
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

            # find quadrent point lies in
            if trans.y_grid[ii,jj] >= massCon.y_grid[nn_ii,nn_jj]:
                i_s = +1
            if trans.x_grid[ii,jj] >= massCon.x_grid[nn_ii,nn_jj]:
                j_s = +1

            # check for missing priority data!
            # NOTE: points are ordered as such:
            #
            # 0: (ii    , jj    )
            # 1: (ii    , jj+j_s)
            # 2: (ii+i_s, jj+j_s)
            # 3: (ii+i_s, jj    )
            #
            # Which, for an upper-right quadrent looks like:
            #
            # 3 ---- 2
            # |      |
            # |      |
            # 0 ---- 1
            #
            # Numbering in other quadrents is the reflection through the axis or axes 
            # with negitive skip values (i_s or j_s).
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
    

