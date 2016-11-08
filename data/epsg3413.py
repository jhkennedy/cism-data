import json
import numpy

from pprint import pprint
from netCDF4 import Dataset

from util.projections import DataGrid

def build_base(f_base, f_epsg, d_meters):
    with open(f_epsg, 'r') as f:
        epsg = json.load(f)

    EPSG = DataGrid()
    EPSG.x = numpy.arange(epsg['ll'][0], epsg['ur'][0], d_meters)
    EPSG.nx = len(EPSG.x)

    EPSG.y = numpy.arange(epsg['ll'][1], epsg['ur'][1], d_meters)
    EPSG.ny = len(EPSG.y)

    nc_base = Dataset(f_base,'w')
    nc_base.createDimension('y',EPSG.ny)
    nc_base.createDimension('x',EPSG.nx)
    
    base = DataGrid()
    base.ny = EPSG.ny
    base.nx = EPSG.nx
    
    base.y = nc_base.createVariable('y', 'f4', 'y')
    base.y[:] = EPSG.y[:]
    #FIXME: attributes

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = EPSG.x[:]
    #FIXME: attributes

    # create some grids for interpolation
    base.make_grid()

    return (nc_base, base)
