from numpy import *
from struct import unpack
from netCDF4 import Dataset
from scipy.interpolate import griddata , interp2d

x0 = -3333500
y0 = -3333500
dx = 1000
dy = 1000
nx = 6667
ny = 6667
file = './originals/1km-res-Bedmap2-DEM/'

print('Creating netCDF file')
nc = Dataset('antarctica_1km.nc', 'w', format='NETCDF4')
x_d  = nc.createDimension('x', nx)
y_d  = nc.createDimension('y', ny)
x_v  = nc.createVariable('x' ,'f4',('x',))
y_v  = nc.createVariable('y' ,'f4',('y',))
#srf_v = nc.createVariable('usrf','f4',('y','x',))
thk_v = nc.createVariable('thk','f4',('y','x',))
topg_v = nc.createVariable('topg','f4',('y','x',))
vx_v = nc.createVariable('vx','f4',('y','x',))
vy_v = nc.createVariable('vy','f4',('y','x',))
verr_v = nc.createVariable('verr','f4',('y','x',))
acab_alb_v = nc.createVariable('acab_alb','f4',('y','x',))
acab_rac_v = nc.createVariable('acab','f4',('y','x',))
artm_alb_v = nc.createVariable('artm_alb','f4',('y','x',))
artm_rac_v = nc.createVariable('artm','f4',('y','x',))
bheatflx_v = nc.createVariable('bheatflx','f4',('y','x',))
#srf_v.missing_value = -9999.
#thk_v.missing_value = -9999.
#topg_v.missing_value = -9999.
acab_rac_v.source = 'racmo, 27km'
x_v[:] = arange(x0,x0+nx*dx,dx)
y_v[:] = arange(y0,y0+ny*dy,dy)

#print('Read bedmap2_surface.flt')
#f = open(file+'bedmap2_bin/bedmap2_surface.flt','rb')
#srf_v[:,:] = asarray( unpack( str(nx*ny)+'f' , f.read() ) ).reshape((ny,nx))[::-1,:]
#f.close()
print('Read bedmap2_thickness.flt')
f = open(file+'bedmap2_bin/bedmap2_thickness.flt','rb')
buf = asarray( unpack( str(nx*ny)+'f' , f.read() ) )
f.close()
buf[buf==-9999.] = 0.
thk_v[:,:] = buf.reshape((ny,nx))[::-1,:]
print('Read bedmap2_bed.flt')
f = open(file+'bedmap2_bin/bedmap2_bed.flt','rb')
topg_v[:,:] = asarray( unpack( str(nx*ny)+'f' , f.read() ) ).reshape((ny,nx))[::-1,:]
f.close()

#################################################################
# VELOCITIES
#################################################################
vel_nx = 6223
vel_ny = 6223
vel_x0 = -2800000
vel_y0 = -2800000
vel_dx = 900
vel_dy = 900
vel_xvals = arange(vel_x0,vel_x0+vel_nx*vel_dx,vel_dx)
vel_yvals = arange(vel_y0,vel_y0+vel_ny*vel_dy,vel_dy)

nc_v = Dataset('originals/Rignot-InSAR-2011/Antarctica_ice_velocity.nc','r')
strd=1
print('Interpolating vx')
f = interp2d( vel_yvals[::strd] , vel_xvals[::strd] , nc_v.variables['vx' ][::strd,::strd] , kind='linear' , copy=False , bounds_error=False, fill_value=0. )
vx_v[:,:]   = f( y_v[:] , x_v[:] )[::-1,:]
print('Interpolating vy')
f = interp2d( vel_yvals[::strd] , vel_xvals[::strd] , nc_v.variables['vy' ][::strd,::strd] , kind='linear' , copy=False , bounds_error=False, fill_value=0. )
vy_v[:,:]   = f( y_v[:] , x_v[:] )[::-1,:]
print('Interpolating err')
f = interp2d( vel_yvals[::strd] , vel_xvals[::strd] , nc_v.variables['err'][::strd,::strd] , kind='linear' , copy=False , bounds_error=False, fill_value=0. )
verr_v[:,:] = f( y_v[:] , x_v[:] )[::-1,:]
nc_v.close()

#################################################################
# BHEATFLX, ACAB, ARTM
#################################################################
alb_nx = 1160
alb_ny = 1120
alb_x0 = -2800000
alb_y0 = -2800000
alb_dx = 5000
alb_dy = 5000
alb_xvals = arange(alb_x0,alb_x0+alb_nx*alb_dx,alb_dx)
alb_yvals = arange(alb_y0,alb_y0+alb_ny*alb_dy,alb_dy)
nc_a = Dataset('originals/ALB-5km-res/ALB_Ant_5km0.nc','r')
strd=1
print('Interpolating acab')
f = interp2d( alb_yvals[::strd] , alb_xvals[::strd] , nc_a.variables['acab'    ][0,::strd,::strd].transpose() , kind='linear' , copy=True  , bounds_error=False, fill_value=0. )
acab_alb_v[:,:] = f( y_v[:] , x_v[:] ).transpose()
print('Interpolating artm')
f = interp2d( alb_yvals[::strd] , alb_xvals[::strd] , nc_a.variables['artm'    ][0,::strd,::strd].transpose() , kind='linear' , copy=False , bounds_error=False, fill_value=0. )
artm_alb_v[:,:]     = f( y_v[:] , x_v[:] ).transpose()
print('Interpolating bheatflx')
f = interp2d( alb_yvals[::strd] , alb_xvals[::strd] , nc_a.variables['bheatflx'][0,::strd,::strd].transpose() , kind='linear' , copy=False , bounds_error=False, fill_value=0. )
bheatflx_v[:,:] = f( y_v[:] , x_v[:] ).transpose()
nc_a.close()

#################################################################
# RACMO: ACAB, ARTM
#################################################################
nc_c = Dataset('coords.nc','r')
xvals = nc_c.variables['xvar'][:]
yvals = nc_c.variables['yvar'][:]
nc_c.close()
print('Interpolating racmo_smb')
nc_r = Dataset('./originals/racmo/1979-2010-average-SMB_AIS.nc','r')
xn,yn = meshgrid( x_v[:] , y_v[:] )
acab_rac_v[:,:] = griddata( ( yvals.flatten() , xvals.flatten() ) , nc_r.variables['smb'][::strd,::strd].flatten() , ( yn.flatten() , xn.flatten() ) , method='linear' , fill_value=0. )
acab_rac_v[:,:] = acab_rac_v[:,:] / 910.
nc_r.close()
nc_r = Dataset('./originals/racmo/1979-2010-average-t2m-AIS.nc','r')
artm_rac_v[:,:] = griddata( ( yvals.flatten() , xvals.flatten() ) , nc_r.variables['t2m'][::strd,::strd].flatten() , ( yn.flatten() , xn.flatten() ) , method='linear' , fill_value=273.15 )
artm_rac_v[:,:] = artm_rac_v[:,:] - 273.15
nc_r.close()


nc.close()
