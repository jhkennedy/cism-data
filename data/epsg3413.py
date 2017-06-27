import json
import numpy as np

from netCDF4 import Dataset

from util import projections


def build_base(f_base, f_epsg, d_meters):
    proj_epsg3413, _ = projections.greenland()
    
    with open(f_epsg, 'r') as f:
        epsg = json.load(f)

    EPSG = projections.DataGrid()
    EPSG.x = np.arange(epsg['ll'][0], epsg['ur'][0], d_meters)
    EPSG.nx = len(EPSG.x)

    EPSG.y = np.arange(epsg['ll'][1], epsg['ur'][1], d_meters)
    EPSG.ny = len(EPSG.y)

    nc_base = Dataset(f_base, 'w', format='NETCDF4')
    nc_base.createDimension('y', EPSG.ny)
    nc_base.createDimension('x', EPSG.nx)

    base = projections.DataGrid()
    base.ny = EPSG.ny
    base.nx = EPSG.nx

    base.y = nc_base.createVariable('y', 'f4', 'y')
    base.y[:] = EPSG.y[:]
    base.y.long_name = 'y-coordinate in projection'
    base.y.standard_name = 'projection_y_coordinate'
    base.y.axis = 'Y'
    base.y.units = 'meters'

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = EPSG.x[:]
    base.x.long_name = 'x-coordinate in projection'
    base.x.standard_name = 'projection_x_coordinate'
    base.x.axis = 'X'
    base.x.units = 'meters'

    base.proj = nc_base.createVariable('epsg_3413', 'b')
    base.proj.grid_mapping_name = 'polar_stereographic'
    base.proj.latitude_of_projection_origin = 90.
    base.proj.straight_vertical_longitude_from_pole = -45.
    base.proj.standard_parallel = 70.
    base.proj.proj_scale_factor = 1.
    base.proj.false_easting = 0.
    base.proj.false_northing = 0.
    base.proj.ellipsoid = 'WGS84'
    base.proj.datum = 'WGS84'
    base.proj.units = 'meters'
    base.proj.proj4_string = projections.proj_string(proj_epsg3413) 

    # create some grids for interpolation
    base.make_grid()

    return (nc_base, base)
