from netCDF4 import Dataset

from util.ncfunc import get_nc_file
from util.ncfunc import copy_atts

lc_mask = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'
nc_mask = get_nc_file(lc_mask, 'r')
y = nc_mask.variables['projection_y_coordinate'][:]
x = nc_mask.variables['projection_x_coordinate'][:]
msk = nc_mask.variables['IceSheetMask'][:,:]

nc_base = Dataset('zurich_mask.nc', 'w')
nc_base.createDimension('x1', 1501)
nc_base.createDimension('y1', 2801)

y1 = nc_base.createVariable('y1', 'f4', ('y1',))
y1.long_name = "y coordinate of projection"
y1.units = "km"
y1.standard_name = "projection_y_coordinate"
y1[:] = y[100:2901]

x1 = nc_base.createVariable('x1', 'f4', ('x1',))
x1.long_name = "x coordinate of projection"
x1.units = "km"
x1.standard_name = "projection_x_coordinate"
x1[:] = x[500:2001]

mask = nc_base.createVariable('mask', 'b', ('y1','x1',))
copy_atts(nc_mask.variables['IceSheetMask'],mask)
mask.units = 'unitless'
mask[:,:] = msk[100:2901,500:2001]

nc_mask.close()
nc_base.close()
