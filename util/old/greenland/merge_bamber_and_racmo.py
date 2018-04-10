
import math
from netCDF4 import Dataset
from numpy import ma,ndarray

def copy_atts(vin,vout) :
  atts = vin.ncattrs()
  for i in range(len(atts)) :
    vout.setncattr(atts[i],vin.getncattr(atts[i]))

nc_b = Dataset( 'Greenland_bedrock_topography_V3.nc' , 'r+' )
nc_o = Dataset( 'Racmo2MeanSMB_1961-1990.nc' , 'r' )

#Read the target grid
x_b  = nc_b.variables['projection_x_coordinate'][:]
y_b  = nc_b.variables['projection_y_coordinate'][:]
nx_b = x_b.shape[0]
ny_b = y_b.shape[0]

#List of variables to merge in
varlist_o = [ 'smb' ]
varlist_b = [ 'acab' ]

vdata_o = ndarray( ( ny_b , nx_b ) )
for v in range(len(varlist_o)) :
  vdata_o[:,:] = 0.
  ncvar_o = nc_o.variables[varlist_o[v]]
  vdata_o[:,:] = ncvar_o[:,::-1].transpose() / 910.
  vdata_o = ma.masked_invalid( vdata_o )
  vdata_o = vdata_o.filled(0.)
  ncvar_b = nc_b.createVariable(varlist_b[v],'f4',('y','x',))
  ncvar_b[:,:] = vdata_o[:,:] 
  copy_atts(ncvar_o,ncvar_b)

nc_b.close()
nc_o.close()
