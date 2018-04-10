from numpy import *
from struct import unpack
from netCDF4 import Dataset

nc500 = Dataset('greenland_vel_mosaic500_bambergrid_reduced.nc', 'r', format='NETCDF4')
nc1km = Dataset('greenland_1km_time.nc', 'r+', format='NETCDF4')
nc2km = Dataset('greenland_2km_time.nc', 'r+', format='NETCDF4')
nc4km = Dataset('greenland_4km_time.nc', 'r+', format='NETCDF4')
nc8km = Dataset('greenland_8km_time.nc', 'r+', format='NETCDF4')

vx500 = nc500.variables['vx'][:,:]
vy500 = nc500.variables['vy'][:,:]
ex500 = nc500.variables['ex'][:,:]
ey500 = nc500.variables['ey'][:,:]

vx_v = nc1km.createVariable('vx','f4',('time','y1','x1',))
vy_v = nc1km.createVariable('vy','f4',('time','y1','x1',))
ex_v = nc1km.createVariable('ex','f4',('time','y1','x1',))
ey_v = nc1km.createVariable('ey','f4',('time','y1','x1',))
vx_v.missing_value = -2.e9
vy_v.missing_value = -2.e9
ex_v.missing_value = -2.e9
ey_v.missing_value = -2.e9
vx_v[0,:,:] = vx500[::2,::2]
vy_v[0,:,:] = vy500[::2,::2]
ex_v[0,:,:] = ex500[::2,::2]
ey_v[0,:,:] = ey500[::2,::2]

vx_v = nc2km.createVariable('vx','f4',('time','y1','x1',))
vy_v = nc2km.createVariable('vy','f4',('time','y1','x1',))
ex_v = nc2km.createVariable('ex','f4',('time','y1','x1',))
ey_v = nc2km.createVariable('ey','f4',('time','y1','x1',))
vx_v.missing_value = -2.e9
vy_v.missing_value = -2.e9
ex_v.missing_value = -2.e9
ey_v.missing_value = -2.e9
vx_v[0,:,:] = vx500[::4,::4]
vy_v[0,:,:] = vy500[::4,::4]
ex_v[0,:,:] = ex500[::4,::4]
ey_v[0,:,:] = ey500[::4,::4]

vx_v = nc4km.createVariable('vx','f4',('time','y1','x1',))
vy_v = nc4km.createVariable('vy','f4',('time','y1','x1',))
ex_v = nc4km.createVariable('ex','f4',('time','y1','x1',))
ey_v = nc4km.createVariable('ey','f4',('time','y1','x1',))
vx_v.missing_value = -2.e9
vy_v.missing_value = -2.e9
ex_v.missing_value = -2.e9
ey_v.missing_value = -2.e9
vx_v[0,:,:] = vx500[::8,::8]
vy_v[0,:,:] = vy500[::8,::8]
ex_v[0,:,:] = ex500[::8,::8]
ey_v[0,:,:] = ey500[::8,::8]

vx_v = nc8km.createVariable('vx','f4',('time','y1','x1',))
vy_v = nc8km.createVariable('vy','f4',('time','y1','x1',))
ex_v = nc8km.createVariable('ex','f4',('time','y1','x1',))
ey_v = nc8km.createVariable('ey','f4',('time','y1','x1',))
vx_v.missing_value = -2.e9
vy_v.missing_value = -2.e9
ex_v.missing_value = -2.e9
ey_v.missing_value = -2.e9
vx_v[0,:,:] = vx500[::16,::16]
vy_v[0,:,:] = vy500[::16,::16]
ex_v[0,:,:] = ex500[::16,::16]
ey_v[0,:,:] = ey500[::16,::16]

nc8km.close()
nc4km.close()
nc2km.close()
nc1km.close()
nc500.close()

