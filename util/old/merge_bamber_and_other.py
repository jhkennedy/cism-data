
from netCDF4 import Dataset
from numpy import ma,ndarray

def copy_atts(vin,vout) :
  atts = vin.ncattrs()
  for i in range(len(atts)) :
    vout.setncattr(atts[i],vin.getncattr(atts[i]))

nc_b = Dataset( 'Greenland_bedrock_topography_V3.nc' , 'r+' )
nc_o = Dataset( 'Greenland1km.nc' , 'r' )

#Read the target grid
x_b  = nc_b.variables['projection_x_coordinate'][:]
y_b  = nc_b.variables['projection_y_coordinate'][:]
nx_b = x_b.shape[0]
ny_b = y_b.shape[0]

#Read the source grid
x_o  = nc_o.variables['x'][:]
y_o  = nc_o.variables['y'][:]
nx_o = x_o.shape[0]
ny_o = y_o.shape[0]

#Find the equal point
xequal = 500
yequal = 100

#List of variables to merge in
#varlist_o = [ 'bheatflx' , 'presartm' , 'smb'  ]
#varlist_b = [ 'bheatflx' , 'artm'     , 'acab' ]
varlist_o = [ 'bheatflx' , 'presartm' ]
varlist_b = [ 'bheatflx' , 'artm'     ]

vdata_o = ndarray( ( ny_b , nx_b ) ) # y-x array, size of projection

# for each variable to merge in
for v in range(len(varlist_o)) :
  # zero everything? in projection array
  vdata_o[:,:] = 0.
  # get variable to merge in...
  ncvar_o = nc_o.variables[varlist_o[v]]
  # starting at [x,y]equal, till end of source grid, set grid variables
  vdata_o[yequal:yequal+ny_o,xequal:xequal+nx_o] = ncvar_o[0,:,:]
  # create a new variable in out file, of type float, with dims y,x and return a variable access class
  ncvar_b = nc_b.createVariable(varlist_b[v],'f4',('y','x',))
  # copy variable merge in data into new variable
  ncvar_b[:,:] = vdata_o[:,:] 
  # inver sign
  if (varlist_o[v] == 'bheatflx') :
    ncvar_b[:,:] = -ncvar_b[:,:]
  # and coppy the global attributes of the variable
  copy_atts(ncvar_o,ncvar_b)

nc_b.close()
nc_o.close()
