import sys
from numpy import *
from struct import unpack
from netCDF4 import Dataset

location = sys.argv[1].rstrip('/') # 0 arg is calling this script, 1 is first input!

print(location)

x0 = -645*1000
y0 = -3370*1000
dx = 500
dy = 500
nx = 3010
ny = 5460
file = location+'/mosaicOffsets'

print('Creating netCDF file')
nc = Dataset(location+'/greenland_vel_mosaic500.nc', 'w', format='NETCDF4')
x_d  = nc.createDimension('x', nx)
y_d  = nc.createDimension('y', ny)
x_v  = nc.createVariable('x' ,'f4',('x',))
y_v  = nc.createVariable('y' ,'f4',('y',))
x_v.units = 'meters'
y_v.units = 'meters'
vx_v = nc.createVariable('vx','f4',('y','x',))
vy_v = nc.createVariable('vy','f4',('y','x',))
ex_v = nc.createVariable('ex','f4',('y','x',))
ey_v = nc.createVariable('ey','f4',('y','x',))
vx_v.units = 'meters/year' 
vy_v.units = 'meters/year' 
ex_v.units = 'meters/year' 
ey_v.units = 'meters/year' 
vx_v.missing_value =-2.e9
vy_v.missing_value =-2.e9
ex_v.missing_value =-2.e9
ey_v.missing_value =-2.e9
vx_v.comments = 'These velocities are true values and not subject to the distance distortions present in a polar stereographic grid.' 
vy_v.comments = 'These velocities are true values and not subject to the distance distortions present in a polar stereographic grid.' 
ex_v.comments = 'These estimates include the statistical uncertainty associated with the phase and speckle tracking error. Formal errors agree reasonably well with errors determined by comparison with GPS data, but underestimate true uncertainty in several ways. They should be used more as an indication of relative quality rather than absolute error.' 
ey_v.comments = 'These estimates include the statistical uncertainty associated with the phase and speckle tracking error. Formal errors agree reasonably well with errors determined by comparison with GPS data, but underestimate true uncertainty in several ways. They should be used more as an indication of relative quality rather than absolute error.' 
vx_v.long_name = 'surface x velocity' 
vy_v.long_name = 'surface y velocity' 
ex_v.long_name = 'surface x velocity error' 
ey_v.long_name = 'surface y velocity error' 
vx_v.standard_name = 'land_ice_surface_x_velocity' 
vy_v.standard_name = 'land_ice_surface_y_velocity' 
vx_v.source = 'Greenland ice sheet surface velocity from InSAR; provided Ian Joughin, April 2015.' 
vy_v.source = 'Greenland ice sheet surface velocity from InSAR; provided Ian Joughin, April 2015.' 
ex_v.source = 'Greenland ice sheet surface velocity from InSAR; provided Ian Joughin, April 2015.' 
ey_v.source = 'Greenland ice sheet surface velocity from InSAR; provided Ian Joughin, April 2015.' 
vx_v.reference = 'Joughin, I., B. Smith, I. Howat, and T. Scambos. 2010. MEaSUREs Greenland Ice Sheet Velocity Map from InSAR Data. [indicate subset used]. Boulder, Colorado USA: National Snow and Ice Data Center.' 
vy_v.reference = 'Joughin, I., B. Smith, I. Howat, and T. Scambos. 2010. MEaSUREs Greenland Ice Sheet Velocity Map from InSAR Data. [indicate subset used]. Boulder, Colorado USA: National Snow and Ice Data Center.' 
ex_v.reference = 'Joughin, I., B. Smith, I. Howat, and T. Scambos. 2010. MEaSUREs Greenland Ice Sheet Velocity Map from InSAR Data. [indicate subset used]. Boulder, Colorado USA: National Snow and Ice Data Center.' 
ey_v.reference = 'Joughin, I., B. Smith, I. Howat, and T. Scambos. 2010. MEaSUREs Greenland Ice Sheet Velocity Map from InSAR Data. [indicate subset used]. Boulder, Colorado USA: National Snow and Ice Data Center.' 

x_v[:] = [i for i in range(x0,x0+nx*dx,dx)]
y_v[:] = [i for i in range(y0,y0+ny*dy,dy)]

print('Reading vx and vy for file')
f = open(file+'.vx','rb')
vx_bin = f.read()
f.close()
f = open(file+'.vy','rb')
vy_bin = f.read()
f.close()
vx = ndarray( shape = ( ny , nx ) , dtype = float )
vy = ndarray( shape = ( ny , nx ) , dtype = float )
vm = ndarray( shape = ( ny , nx ) , dtype = float )
for j in range(ny) :
  if (j%200 == 0) :
    print("Percentage complete: ",float(j)/ny*100.)
  for i in range(nx) :
    ind = j*nx+i
    beg =  ind   *4
    end = (ind+1)*4
    vx[j,i] = float( unpack('>f',vx_bin[beg:end])[0] )
    vy[j,i] = float( unpack('>f',vy_bin[beg:end])[0] )

print('writing vx, vy, and vm to file')
vx_v[:,:] = vx
vy_v[:,:] = vy
print(vx.max())
print(vy.max())

print('Reading ex and ey')
f = open(file+'.ex','rb')
ex_bin = f.read()
f.close()
f = open(file+'.ey','rb')
ey_bin = f.read()
f.close()
ex = ndarray( shape = ( ny , nx ) , dtype = float )
ey = ndarray( shape = ( ny , nx ) , dtype = float )
em = ndarray( shape = ( ny , nx ) , dtype = float )
for j in range(ny) :
  if (j%200 == 0) :
    print("Percentage complete: ",float(j)/ny*100.)
  for i in range(nx) :
    ind = j*nx+i
    beg =  ind   *4
    end = (ind+1)*4
    ex[j,i] = float( unpack('>f',ex_bin[beg:end])[0] )
    ey[j,i] = float( unpack('>f',ey_bin[beg:end])[0] )
    if ( ex[j,i] > 1.e9 or ey[j,i] > 1.e9 ) :
      ex[j,i] = -2e9
      ey[j,i] = -2e9

print('writing ex, ey, and em to file')
ex_v[:,:] = ex
ey_v[:,:] = ey
print(ex.max())
print(ey.max())

nc.close()

