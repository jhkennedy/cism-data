import os
import scipy
import pyproj
import numpy as np
from scipy.spatial import cKDTree

from util.ncfunc import get_nc_file
from util.projections import DataGrid

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
massCon.make_grid()

massCon.yx = np.ndarray( (len( massCon.y_grid.ravel() ),2) )
massCon.yx[:,0] = massCon.y_grid.ravel()
massCon.yx[:,1] = massCon.x_grid.ravel()

# ckd trees of data grids
bamber.tree = cKDTree(bamber.yx)
massCon.tree = cKDTree(massCon.yx)

# get transformation grids
proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')
proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m')

b_trans = DataGrid()
b_trans.ny = bamber.ny
b_trans.nx = bamber.nx

b_trans.x_grid, b_trans.y_grid = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, bamber.x_grid.flatten(), bamber.y_grid.flatten())
b_trans.y_grid = b_trans.y_grid.reshape((b_trans.ny,b_trans.nx))
b_trans.x_grid = b_trans.x_grid.reshape((b_trans.ny,b_trans.nx))

b_trans.yx = np.ndarray( (len( b_trans.y_grid.ravel() ),2) )
b_trans.yx[:,0] = b_trans.y_grid.ravel()
b_trans.yx[:,1] = b_trans.x_grid.ravel()

m_trans = DataGrid()
m_trans.ny = massCon.ny
m_trans.nx = massCon.nx

m_trans.x_grid, m_trans.y_grid = pyproj.transform(proj_epsg3413, proj_eigen_gl04c, massCon.x_grid.flatten(), massCon.y_grid.flatten())
m_trans.y_grid = m_trans.y_grid.reshape((m_trans.ny,m_trans.nx))
m_trans.x_grid = m_trans.x_grid.reshape((m_trans.ny,m_trans.nx))

m_trans.yx = np.ndarray( (len( m_trans.y_grid.ravel() ),2) )
m_trans.yx[:,0] = m_trans.y_grid.ravel()
m_trans.yx[:,1] = m_trans.x_grid.ravel()

# find transform grids nearest neighbor
b_trans.qd, b_trans.qi = massCon.tree.query(b_trans.yx, k=1)
m_trans.qd, m_trans.qi = bamber.tree.query(m_trans.yx, k=1)

# build data array.
new_data = np.ma.array( np.zeros( (bamber.ny,bamber.nx) ), mask=np.zeros( (bamber.ny,bamber.nx) ) )


