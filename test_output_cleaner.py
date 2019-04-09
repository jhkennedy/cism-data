#!/usr/bin/env python
# encoding: utf-8

import scipy
import pyproj

from netCDF4 import Dataset

FILE = 'test_output.nc'


with Dataset(FILE, 'r+') as nc:
    nc.variables['y'].units = 'meters'
    nc.variables['x'].units = 'meters'

    # nc.variables['z'].grid_mapping = 'crs'  # already an attribute
    nc.variables['z'].coordinates = 'lon lat'

    # String generated from `$ gdalsrsinfo test_output.prj`
    # test_output.prj created from nc.variables['crs'].spatial_ref
    nc.variables['crs'].proj4_str = '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 ' \
                                    '+x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 ' \
                                    '+units=m +no_defs'
    # nc.variables['crs'].grid_mapping_name = 'albers_coords'  # already an attribute
    nc.variables['crs'].latitude_of_projection_origin = 50.
    nc.variables['crs'].straight_vertical_longitude_from_pole = -154.
    nc.variables['crs'].standard_parallel = (55., 65.)
    nc.variables['crs'].proj_scale_factor = 1.
    nc.variables['crs'].false_easting = 0.
    nc.variables['crs'].false_northing = 0.
    nc.variables['crs'].ellipsoid = 'GSR80'
    nc.variables['crs'].datum = 'NAD83'
    nc.variables['crs'].towgs84 = '0,0,0,0,0,0,0'
    # nc.variables['crs'].units = 'm'  # already an attribute

    y = nc.variables['y'][:]
    x = nc.variables['x'][:]
    nc.variables['y'][:] = y.round()
    nc.variables['x'][:] = x.round()

    proj = pyproj.Proj(nc.variables['crs'].proj4_str)
    y_grid, x_grid = scipy.meshgrid(y.round(), x.round(), indexing='ij')

    nc_lon_grid = nc.createVariable('lon', 'f4', ('y', 'x',))
    nc_lon_grid.long_name = 'grid center longitude'
    nc_lon_grid.standard_name = 'longitude'
    nc_lon_grid.units = 'degrees_east'
    nc_lon_grid.grid_mapping = nc.variables['z'].grid_mapping

    nc_lat_grid = nc.createVariable('lat', 'f4', ('y', 'x',))
    nc_lat_grid.long_name = 'grid center latitude'
    nc_lat_grid.standard_name = 'latitude'
    nc_lat_grid.units = 'degrees_north'
    nc_lat_grid.grid_mapping = nc.variables['z'].grid_mapping

    lon_grid, lat_grid = proj(x_grid.ravel(), y_grid.ravel(), inverse=True)
    lon_grid.shape = x_grid.shape
    lat_grid.shape = x_grid.shape

    nc_lon_grid[:, :] = lon_grid
    nc_lat_grid[:, :] = lat_grid


