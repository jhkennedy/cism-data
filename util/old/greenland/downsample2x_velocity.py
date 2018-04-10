
from numpy import *
from struct import unpack
from netCDF4 import Dataset

def downsample( a , f ) :
  if (f == 1) :
    return a
  nx = a.shape[0] / f
  ny = a.shape[1] / f
  a1 = a[:a.shape[0],:ny]
  for j in range(ny) :
    a1[:,j] = a[:,j*f:(j+1)*f].mean(1)
  a2 = a1[:nx,:ny]
  for i in range(nx) :
    a2[i,:] = a1[i*f:(i+1)*f,:].mean(0)
  return a2

nc500 = Dataset('greenland_vel_mosaic500_bambergrid.nc', 'r', format='NETCDF4')
nc1km = Dataset('greenland_vel_mosaic1km_bambergrid.nc', 'w', format='NETCDF4')

x500 = nc500.variables['x'][:]
y500 = nc500.variables['y'][:]
vx500 = nc500.variables['vx'][:]
vy500 = nc500.variables['vy'][:]
ex500 = nc500.variables['ex'][:]
ey500 = nc500.variables['ey'][:]

nx500 = x500.shape[0]
ny500 = y500.shape[0]

nx1km = nx500/2
ny1km = ny500/2

x1km  = ndarray((nx1km))
y1km  = ndarray((ny1km))
vx1km = ndarray((ny1km,nx1km))
vy1km = ndarray((ny1km,nx1km))
ex1km = ndarray((ny1km,nx1km))
ey1km = ndarray((ny1km,nx1km))

for i in range(nx1km) :
  x1km[i] = ( x500[2*i] + x500[2*i+1] ) / 2.

for j in range(ny1km) :
  y1km[j] = ( y500[2*j] + y500[2*j+1] ) / 2.

print(y500)
print(y1km)

vx1km = downsample( vx500 , 2 )
vy1km = downsample( vy500 , 2 )
ex1km = downsample( ex500 , 2 )
ey1km = downsample( ey500 , 2 )

x_d  = nc1km.createDimension('x', nx1km)
y_d  = nc1km.createDimension('y', ny1km)
x_v  = nc1km.createVariable('x' ,'f4',('x',))
y_v  = nc1km.createVariable('y' ,'f4',('y',))
vx_v = nc1km.createVariable('vx','f4',('y','x',))
vy_v = nc1km.createVariable('vy','f4',('y','x',))
ex_v = nc1km.createVariable('ex','f4',('y','x',))
ey_v = nc1km.createVariable('ey','f4',('y','x',))
vx_v.missing_value = -2.e9
vy_v.missing_value = -2.e9
ex_v.missing_value = -2.e9
ey_v.missing_value = -2.e9
x_v[:] = x1km
y_v[:] = y1km
vx_v[:,:] = vx1km
vy_v[:,:] = vy1km
ex_v[:,:] = ex1km
ey_v[:,:] = ey1km

nc1km.close()
nc500.close()

