
from netCDF4 import Dataset
from numpy import ma,ndarray

def copy_atts(vin,vout) :
  atts = vin.ncattrs()
  for i in range(len(atts)) :
    vout.setncattr(atts[i],vin.getncattr(atts[i]))

nc_1 = Dataset( 'greenland_1km.nc' , 'r+' )
nc_2 = Dataset( 'greenland_2km.nc' , 'w'  )
nc_4 = Dataset( 'greenland_4km.nc' , 'w'  )
nc_8 = Dataset( 'greenland_8km.nc' , 'w'  )

#Read the target grid
x_1  = nc_1.variables['x'][:]
y_1  = nc_1.variables['y'][:]
nx_1 = x_1.shape[0]
ny_1 = y_1.shape[0]

x_2 = x_1[::2]
y_2 = y_1[::2]
x_4 = x_1[::4]
y_4 = y_1[::4]
x_8 = x_1[::8]
y_8 = y_1[::8]

dx2 = nc_2.createDimension('x', len(x_2))
dy2 = nc_2.createDimension('y', len(y_2))
dx4 = nc_4.createDimension('x', len(x_4))
dy4 = nc_4.createDimension('y', len(y_4))
dx8 = nc_8.createDimension('x', len(x_8))
dy8 = nc_8.createDimension('y', len(y_8))

vx2 = nc_2.createVariable('x','f4',('x',))
vy2 = nc_2.createVariable('y','f4',('y',))
vx4 = nc_4.createVariable('x','f4',('x',))
vy4 = nc_4.createVariable('y','f4',('y',))
vx8 = nc_8.createVariable('x','f4',('x',))
vy8 = nc_8.createVariable('y','f4',('y',))

vx2[:] = x_2
vy2[:] = y_2
vx4[:] = x_4
vy4[:] = y_4
vx8[:] = x_8
vy8[:] = y_8

#List of variables to merge in
varlist = [ 'acab' , 'artm' , 'bheatflx' , 'thk' , 'topg' , 'topgerr' , 'usrf' , 'usrfRMSE'  ]

for v in range(len(varlist)) :
  vx2 = nc_2.createVariable(varlist[v],'f4',('y','x',))
  vx4 = nc_4.createVariable(varlist[v],'f4',('y','x',))
  vx8 = nc_8.createVariable(varlist[v],'f4',('y','x',))
  vx2[:,:] = nc_1.variables[varlist[v]][::2,::2]
  vx4[:,:] = nc_1.variables[varlist[v]][::4,::4]
  vx8[:,:] = nc_1.variables[varlist[v]][::8,::8]
  copy_atts( nc_1.variables[varlist[v]] , vx2 )
  copy_atts( nc_1.variables[varlist[v]] , vx4 )
  copy_atts( nc_1.variables[varlist[v]] , vx8 )

nc_1.close()
nc_2.close()
nc_4.close()
nc_8.close()

