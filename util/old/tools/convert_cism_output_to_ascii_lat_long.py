# Script to take CISM output and convert to ascii files of the format:
# latitude, longitude, thickness
# This involves reprojecting from polar stereographic to lat long
# and then building 1-d arrays to save as text files for each time slice.

import netCDF4
import numpy as np
import pyproj
import sys


####################### SET THIS STUFF ###################

infile = './processed.out.nc'  # MH
#infile = '/Users/sprice/work/data/piscees/GreenlandFluxSMBforcing/greenland-flux-forcing-simulation/post_process/processed.out.nc' # SP

output_projection = 'LL'   # LL or GIMP or Bamber

outputType = 'thk'  # can be 'thk' (for GRACE) or 'usrf' (for ICESat)

cells_with_smb_only = 1  # 1=only keep cells that have ice AND RACMO SMB; 0=keep all cells with ice
                         # SFP: for newest runs, we don't need to do this step, since SMB exists everywhere
                         # thickness exists
maskfile = './gis4km.smbmask.nc'  # MH
#maskfile = '/Users/mhoffman/documents/greenland_flux_gates/setup_4km_SMB_forcing/gis4km.smbmask.nc'  # MH
##########################################



# Get stuff from main file
f = netCDF4.Dataset(infile,'r')
x1 = f.variables['x1'][:]
y1 = f.variables['y1'][:]
thk = f.variables['thk']
topg = f.variables['topg']

###usrf = thk[:] + topg[:]  # doing this is probably a bad idea because the whole array must be in memory

nx = x1.size
ny = y1.size

nt = len(f.dimensions['time'])

ncells = nx * ny


# Get smbmask, if needed
if cells_with_smb_only == 1:
  fmask = netCDF4.Dataset(maskfile,'r')
  smbmask = fmask.variables['smbmask']


# =====================
# setup proj object
# =====================

#======== Notes about pyproj =========
# The documentation for the transform function says:"In addition to converting between cartographic and geographic projection coordinates, this function can take care of datum shifts (which cannot be done using the __call__ method of the Proj instances). "
# Therefore I am switching from the call method of the Proj instance (used previously) to the transform function so that we can include datum shifts.
# Since this becomes an x,y,z transform instead of just x,y, we can't do it once prior to the time step loop - we must now call the function for every x,y position at each time, since the z will change on each time.  This will make the script a lot slower.
# ===================================

# projection = "Polar-stereographic with standard parallel 71 deg N and central longitude 321 deg E" ;
#'+proj=stere +lat_ts=Latitude at natural origin 
#              +lat_0=90
#              +lon_0=Longitude at natural origin
#	      +k_0=Scale factor at natural origin (normally 1.0)
#              +x_0=False Easting
#              +y_0=False Northing'

#shift_x = -800000           #added to x1 before map conversion
#shift_y = -3400000          #added to y1 before map conversion

# the larger ice2sea grid has additional shifts of 500,000 in x and 100,000 in y
# the offsets should be negative when used as false easting/northing.

# Below line is for datasets used with the 5 km res. ice2sea runs
#projCISM = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=1300000.0 +y_0=3500000.0 +geoidgrids=/Users/mhoffman/documents/greenland_flux_gates/post_process/egm08_25.gtx')  # includes datum shift.  The .gtx file needs to be downloaded and the path adjusted for where you save it.

# Below line is for the newer, 4 km PISCEES-generated datasets 
#projBamber_nogeoid = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0') # OLD VERSION WITHOUT GEOID SPECIFICATION 


# ======== DEFINE PROJECTIONS =============
# CISM's projection is as follows, with the vertical datum as EIGEN-GL04C geoid. 
# datum is actually EIGEN-GL04C but that is not an option in Proj.  Therefore using EGM08 which should be within ~1m everywhere (and 10-20 cm in most places)
# NOTE!!!!!!  egm08_25.gtx can be downloaded from:  http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx  and the path in the projection specification line should point to it!
projBamber = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids=./egm08_25.gtx')

# Standard Lat/Long 
projLL   = pyproj.Proj(proj='latlong', datum='WGS84')

# GIMP projection: This is also polar stereographic but with different standard parallel and using the WGS84 ellipsoid.
projGIMP = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=315.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
# ===================================



# For testing, switch all the lines like np.arange(0,nx,1) to something like np.arange(0,nx,10).


# determine lat long values for each point in the grid
# ncells is the max possible size if ice is everywhere, but the actual size will end up being less.
xIn = np.zeros((ncells,1))
yIn = np.zeros((ncells,1))
zIn = np.zeros((ncells,1))
xOut = np.zeros((ncells,1))
yOut = np.zeros((ncells,1))
zOut = np.zeros((ncells,1))
# create new data arrays that are collapses to 1-d array to match lat long arrays
usrfColumn = np.zeros((ncells,1))
thkColumn = np.zeros((ncells,1))


if output_projection == 'LL':
   projOut = projLL
elif output_projection == 'GIMP':
   projOut = projGIMP
elif output_projection == 'Bamber':
   projOut = projBamber
else:
   sys.exit('Error: unknown projection specified!')

#for t in range(15):  # range(nt)
#  time = 2000 + t
for t in range(nt):         #SP: realign start time; allow for half year increments 
#  time = 1991.0 + t/2.0
  time = 1991.0 + t
  print '==== Doing year', time
  print 'each . printed is 100 records completed:'

  ind = 0
  for i in np.arange(0,nx,1):
    for j in np.arange(0,ny,1):
      thkHere = thk[t,j,i]
      if thkHere > 0.0:  # only bother if there is ice
       if (cells_with_smb_only == 1 and smbmask[0,j,i] == 1) or cells_with_smb_only == 0:  # if only wanting cells with smb, then make sure we have smb!
        if ind%100 == 0:
           sys.stdout.write('.'); sys.stdout.flush()  # print a dot for every 100 records
        # convert x,y,z to Lat,Long,Z
        zHere = thkHere + topg[t,j,i]  # This is usrf but we might not have that.  Assuming grounded ice!!!!  
        thkColumn[ind] = thkHere   # This is only needed for GRACE

        xIn[ind] = x1[i]
        yIn[ind] = y1[j]
        zIn[ind] = zHere

        if output_projection == 'Bamber':  # no projection needed
           xOut[ind] = x1[i]; yOut[ind] = y1[j]; zOut[ind] = zHere
        else:
           xOut[ind], yOut[ind], zOut[ind] = pyproj.transform(projBamber, projOut, x1[i], y1[j], z=zHere)
        ###lons[ind], lats[ind] = projBamber_nogeoid(x1[i], y1[i], inverse=True); llZ[ind] = zHere    # OLD METHOD that does not include datum shift
        ind = ind + 1

  print '\nCells used (those with ice):', ind+1
  print 'Cells skipped (those without ice):', ncells - (ind+1)


  # write file for this time
  # Below, the [:ind] notation will write out everything through the last active index
  # ================ ICESat =====================
  if outputType == 'usrf':
      fname="cism_usrf_yr_%f.txt"%(time)        #SP: allow for decimal yrs in output format
      ## write out lat, lon, elev (for NASA altimetry comparisons)
      if output_projection == 'LL':
          # NASA (Jack Saba) currently are requesting format of: Lat, Long, wgs84elev, cism_x, cism_y, cism_z
          np.savetxt(fname,np.concatenate((yOut[:ind], xOut[:ind], zOut[:ind], xIn[:ind], yIn[:ind], zIn[:ind]),axis=1), fmt='%.8f, %.8f, %.2f, %.2f, %.2f, %.2f', delimiter=', ', header='latitude, longitude, surf_elev_wgs84, cism_x, cism_y, surf_elev_geoid')
      else:
          np.savetxt(fname,np.concatenate((xOut[:ind], yOut[:ind], zOut[:ind]),axis=1), fmt='%12.8f', delimiter=', ', header='x, y, surf_elev_m')

  # ================ GRACE =====================
  elif outputType == 'thk':
      fname="cism_thk_yr_%f.txt"%(time)        #SP: allow for decimal yrs in output format
      ## write out lat, lon, thk (for GRACE gravity comparisons)
      if output_projection == 'LL':
          np.savetxt(fname,np.concatenate((yOut[:ind], xOut[:ind], thkColumn[:ind]),axis=1), fmt='%.8f, %.8f, %.2f', delimiter=', ', header='latitude, longitude, thickness_m')
      else:
          np.savetxt(fname,np.concatenate((xOut[:ind], yOut[:ind], thkColumn[:ind]),axis=1), fmt='%12.8f', delimiter=', ', header='x, y, thickness_m')




## write out lat, lon, thick (for GRACE gravity omparisons)
#  np.savetxt(fname,np.concatenate((lats[keepCells], lons[keepCells], thkColumn[keepCells]),axis=1), fmt='%12.8f', delimiter=', ', header='latitude, longitude, thickness_m')

## write output w/ polar stereo x,y rather than lat, lon
#  np.savetxt(fname,np.concatenate((lats[keepCells], lons[keepCells], x[keepCells], y[keepCells], dataColumn[keepCells]),axis=1), fmt='%12.8f', delimiter=', ', header='lat, lon, x, y, surf_elev_m')

##Optional - plot out the lat long positions to make sure they look ok
import matplotlib.pyplot as plt
fig = plt.figure(1)
plt.scatter(xOut[:ind], yOut[:ind], 12, thkColumn[:ind], edgecolors='none')
plt.colorbar()
plt.show()


