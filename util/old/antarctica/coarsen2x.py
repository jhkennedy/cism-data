
from netCDF4 import Dataset
from numpy import ma,ndarray,meshgrid
from scipy.interpolate import griddata , interp2d

def copy_atts(vin,vout) :
  atts = vin.ncattrs()
  for i in range(len(atts)) :
    vout.setncattr(atts[i],vin.getncattr(atts[i]))

nc_h = Dataset( 'antarctica_500m.nc' , 'r' )
nc_1 = Dataset( 'antarctica_1km.nc'  , 'w' )
nc_2 = Dataset( 'antarctica_2km.nc'  , 'w'  )
nc_4 = Dataset( 'antarctica_4km.nc'  , 'w'  )
nc_8 = Dataset( 'antarctica_8km.nc'  , 'w'  )

#Read the target grid
x_h  = nc_h.variables['x'][:]
y_h  = nc_h.variables['y'][:]
nx_h = x_h.shape[0]
ny_h = y_h.shape[0]

x_1 = x_h[::2]
y_1 = y_h[::2]
x_2 = x_h[::4]
y_2 = y_h[::4]
x_4 = x_h[::8]
y_4 = y_h[::8]
x_8 = x_h[::16]
y_8 = y_h[::16]

nx_1 = len(x_1)
ny_1 = len(y_1)
nx_2 = len(x_2)
ny_2 = len(y_2)
nx_4 = len(x_4)
ny_4 = len(y_4)
nx_8 = len(x_8)
ny_8 = len(y_8)

dx1 = nc_1.createDimension('x', len(x_1))
dy1 = nc_1.createDimension('y', len(y_1))
dx2 = nc_2.createDimension('x', len(x_2))
dy2 = nc_2.createDimension('y', len(y_2))
dx4 = nc_4.createDimension('x', len(x_4))
dy4 = nc_4.createDimension('y', len(y_4))
dx8 = nc_8.createDimension('x', len(x_8))
dy8 = nc_8.createDimension('y', len(y_8))

vx1 = nc_1.createVariable('x','f4',('x',))
vy1 = nc_1.createVariable('y','f4',('y',))
vx2 = nc_2.createVariable('x','f4',('x',))
vy2 = nc_2.createVariable('y','f4',('y',))
vx4 = nc_4.createVariable('x','f4',('x',))
vy4 = nc_4.createVariable('y','f4',('y',))
vx8 = nc_8.createVariable('x','f4',('x',))
vy8 = nc_8.createVariable('y','f4',('y',))

vx1[:] = x_1
vy1[:] = y_1
vx2[:] = x_2
vy2[:] = y_2
vx4[:] = x_4
vy4[:] = y_4
vx8[:] = x_8
vy8[:] = y_8

#List of variables to merge in
varlist = [ 'acab' , 'acab_alb' , 'artm' , 'artm_alb' , 'bheatflx' , 'thk' , 'topg' , 'vx' , 'vy' , 'verr' ]

strd = 1

for v in range(len(varlist)) :
  print(varlist[v])
  vx1 = nc_1.createVariable(varlist[v],'f4',('y','x',))
  vx2 = nc_2.createVariable(varlist[v],'f4',('y','x',))
  vx4 = nc_4.createVariable(varlist[v],'f4',('y','x',))
  vx8 = nc_8.createVariable(varlist[v],'f4',('y','x',))
  print('500m')
  f = interp2d( y_h[::strd] , x_h[::strd] , nc_h.variables[varlist[v]][::strd,::strd] , kind='linear' , copy=False , bounds_error=False, fill_value=0. )
  print('1km')
  vx1[:,:] = f( y_1 , x_1 )
  print('2km')
  vx2[:,:] = f( y_2 , x_2 )
  print('4km')
  vx4[:,:] = f( y_4 , x_4 )
  print('8km')
  vx8[:,:] = f( y_8 , x_8 )
  copy_atts( nc_h.variables[varlist[v]] , vx1 )
  copy_atts( nc_h.variables[varlist[v]] , vx2 )
  copy_atts( nc_h.variables[varlist[v]] , vx4 )
  copy_atts( nc_h.variables[varlist[v]] , vx8 )

nc_h.close()
nc_1.close()
nc_2.close()
nc_4.close()
nc_8.close()

