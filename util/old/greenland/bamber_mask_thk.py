

from netCDF4 import Dataset
from numpy import ma,ndarray

def copy_atts(vin,vout) :
  atts = vin.ncattrs()
  for i in range(len(atts)) :
    vout.setncattr(atts[i],vin.getncattr(atts[i]))

nc_b = Dataset( 'Greenland_bedrock_topography_V3.nc'                 , 'r+' )
nc_m = Dataset( 'ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc' , 'r'  )

#Read the target grid
x_b  = nc_b.variables['x'][:]
y_b  = nc_b.variables['y'][:]
nx_b = x_b.shape[0]
ny_b = y_b.shape[0]

v_thk = nc_b.variables['thk']
topg  = nc_b.variables['topg'][:,:].filled(0.)
nc_b.variables['topg'][:,:] = topg
mask  = nc_m.variables['IceSheetMask'][:,:]
thk = v_thk[:,:].filled(0.)
thk = thk * mask
v_thk[:,:] = thk
v_usrf = nc_b.createVariable('usrf','f4',('y','x',))
v_usrf[:,:] = (topg + thk)

nc_b.close()
nc_m.close()

