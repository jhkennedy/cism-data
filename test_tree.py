
import numpy as np
import scipy

from util.ncfunc import get_nc_file
from util.projections import DataGrid

lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'


nc_bamber = get_nc_file(lc_bamber,'r')
nc_massCon = get_nc_file(lc_massCon,'r')

bamber  = DataGrid()
massCon = DataGrid()

bamber.y = nc_bamber.variables['projection_y_coordinate']
bamber.x = nc_bamber.variables['projection_x_coordinate']
bamber.ny = bamber.y[:].shape[0]
bamber.nx = bamber.x[:].shape[0]
bamber.make_grid()

massCon.y = nc_massCon.variables['y']
massCon.x = nc_massCon.variables['x']
massCon.ny = massCon.y[:].shape[0]
massCon.nx = massCon.x[:].shape[0]
massCon.make_grid()

bamber.r_grid = ( bamber.y_grid.ravel(), bamber.x_grid.ravel() )
