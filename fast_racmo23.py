#!/usr/bin/env python2
# encoding: utf-8

"""
Build the RAMCO 2.3 Greenland monthly data on the Bamber grid quick and dirty like
"""

import json
import scipy.interpolate
import numpy as np

from netCDF4 import Dataset

from util import ncfunc
from util import projections


nc_racmo = ncfunc.get_nc_file('data/RACMO2.3/monthly/smb.1958-2013.BN_1958_2013.MM.nc', 'r')
racmo = projections.DataGrid()
racmo.time = nc_racmo.variables['time']
racmo.lat = nc_racmo.variables['lat']
racmo.lon = nc_racmo.variables['lon']
racmo.lat_grid = nc_racmo.variables['LAT']  # Dims: lat, lon
racmo.lon_grid = nc_racmo.variables['LON']  # Dims: lat, lon
racmo.smb = nc_racmo.variables['smb']  # Dims: time, lat, lon


with open('data/BamberDEM/Bambergrid_shrunk.json', 'r') as f:
    bamber = json.load(f)

km_reso = 8
nc_new = Dataset('smb.1958-2013.BN.MM.bamber-{}km.nc'.format(km_reso), 'w')
new = projections.DataGrid()
y_1km = np.arange(bamber['ll'][1], bamber['ur'][1]+1000., 1000.)
x_1km = np.arange(bamber['ll'][0], bamber['ur'][0]+1000., 1000.)
new.y = y_1km[::km_reso]
new.x = x_1km[::km_reso]
new.ny = new.y.shape[0]
new.nx = new.x.shape[0]
new.make_grid()

nc_new.createDimension('time', None)
nc_new.createDimension('y1', new.ny)
nc_new.createDimension('x1', new.nx)

new.time = nc_new.createVariable('time', 'f4', ('time',))
new.time[:] = racmo.time[:]
new.time.long_name = 'time'
new.time.standard_name = 'time'
new.time.units = racmo.time.units
new.time.calendar = '365_day'
# FIXME new.time.comments = FIXME

new.y1 = nc_new.createVariable('y1', 'f4', ('y1',))
new.y1[:] = new.y[:]
new.y1.long_name = 'y-coordinate of projection'
new.y1.standard_name = 'projection_y_coordinate'
# new.y1.axis = 'Y'
new.y1.units = 'm'

new.x1 = nc_new.createVariable('x1', 'f4', ('x1',))
new.x1[:] = new.x[:]
new.x1.long_name = 'x-coordinate of projection'
new.x1.standard_name = 'projection_x_coordinate'
# new.x1.axis = 'X'
new.x1.units = 'm'

proj_name = 'BamberGrid'
_, eigenGl04c = projections.greenland()
nc_new, new = projections.grid_center_latlons(nc_new, new, eigenGl04c, proj_name)

nc_new.createVariable(proj_name, 'c')
nc_new.variables[proj_name].grid_mapping_name = "polar_stereographic"
nc_new.variables[proj_name].straight_vertical_longitude_from_pole = 321.
nc_new.variables[proj_name].latitude_of_projection_origin = 90.
nc_new.variables[proj_name].standard_parallel = 71.
nc_new.variables[proj_name].false_easting = 0.
nc_new.variables[proj_name].false_northing = 0.
nc_new.variables[proj_name].scale_factor_at_projection_origin = 1.
nc_new.variables[proj_name].units = "m"
nc_new.variables[proj_name].proj4_string = '+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0'

racmo.xs, racmo.ys = eigenGl04c(racmo.lon_grid[:, :].ravel(), racmo.lat_grid[:, :].ravel())  # pyproj is gross and does (lon, lat) for in and out. Barf.

new.smb = nc_new.createVariable('smb', 'f4', ('time', 'y1', 'x1',))
new.smb.long_name = "water equivalent surface mass balance"
new.smb.standard_name = "land_ice_lwe_surface_specific_mass_balance"
new.smb.units = "mm year-1"
new.smb.grid_mapping = proj_name
new.smb.coordinates = "lon lat"
new.smb.source = "J. T. M. Lenaerts"
new.smb.reference = "Noel, B., van de Berg, W. J., Machguth, H., Lhermitte, S., Howat, I., Fettweis, X., and van den Broeke, M. R.: a daily, 1 km resolution data set of downscaled Greenland ice sheet surface mass balance (1958--2015), The Cryosphere, 10, 2361-2377, doi:10.5194/tc-10-2361-2016, 2016."
new.smb.comments = "2D linear interpolation of 11 km dataset."

for ii in range(racmo.time[:].shape[0]):
    new.smb[ii, :, :] = scipy.interpolate.griddata(zip(racmo.ys, racmo.xs), racmo.smb[ii, :, :].ravel(),
                                                   (new.y_grid, new.x_grid), method='linear')


nc_new.close()
