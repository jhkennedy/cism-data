
from netCDF4 import Dataset
from numpy import ma,ndarray

def copy_atts(vin,vout) :
  atts = vin.ncattrs()
  for i in range(len(atts)) :
    vout.setncattr(atts[i],vin.getncattr(atts[i]))


def add_time_dim(fname) :
  nc_o = Dataset( fname+'.nc'      , 'r' )
  nc_n = Dataset( fname+'_time.nc' , 'w' )
  
  yo = nc_o.variables['y']
  xo = nc_o.variables['x']
  
  tdn = nc_n.createDimension('time',None   )
  ydn = nc_n.createDimension('y1',len(yo))
  xdn = nc_n.createDimension('x1',len(xo))
  
  tn = nc_n.createVariable('time','f4',('time',))
  yn = nc_n.createVariable('y1','f4',('y1',))
  xn = nc_n.createVariable('x1','f4',('x1',))
  
  xn[:] = xo[:]
  yn[:] = yo[:]
  tn[0] = 0.
  
  varlist = [ 'acab' , 'acab_alb' , 'artm' , 'artm_alb' , 'bheatflx' , 'thk' , 'topg' , 'vx' , 'vy' , 'verr' ]
  
  for v in range(len(varlist)) :
    ncvar_n = nc_n.createVariable(varlist[v],'f4',('time','y1','x1',))
    ncvar_o = nc_o.variables[varlist[v]]
    ncvar_n[0,:,:] = ncvar_o[:,:]
    copy_atts(ncvar_o,ncvar_n)
  
  nc_n.close()
  nc_o.close()

add_time_dim('antarctica_500m')
add_time_dim('antarctica_1km')
add_time_dim('antarctica_2km')
add_time_dim('antarctica_4km')
add_time_dim('antarctica_8km')


