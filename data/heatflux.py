"""
data.heatflux : Antarctic basal heat flux import module.

This module provides ...

Functions list:
    * ...

Notes
-----

The data's projection is specified as "polar steographic", but no information
given to indicate the standard parallel (true scale latitude). Therefore, I am
ASSUMING the projection information is such:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = -71 degrees
    * Latitude of projection origin = -90 degrees
    * Central meridian = 0 degrees
    * false eastings = 0
    * flase northings = 0


However, the standard parallel could be -90 degrees (pole), -81d06m52.3s (ArcGIS
default), -70 degrees (EPSG:3412; often used by NSIDC) -- -71 degrees (EPSG:3031)
does seem to be the most common for glaciology though.

The data was provided as an ascii text file, with a list of x,y point in the
projection (whatever that is), and a z value indicating the heat flux in mW/m^2.

References
----------
Martos, Yasmina M (2017): Antarctic geothermal heat flux distribution and
estimated Curie Depths, links to gridded files. PANGAEA,
https://doi.org/10.1594/PANGAEA.882503.

Supplement to:

Martos, Yasmina M; Catalan, Manuel; Jordan, Tom A; Golynsky, Alexander V;
Golynsky, Dmitry A; Eagles, Graeme; Vaughan, David G (2017): Heat flux
distribution of Antarctica unveiled. Geophysical Research Letters, 44(22),
11417-11426, https://doi.org/10.1002/2017GL075609.

License:

Dataset is licensed under the Creative Commons Attribution 3.0 Unported license
http://creativecommons.org/licenses/by/3.0/
"""

import os
import scipy.interpolate

import numpy as np



def bheatflx(args, nc_base, base):
    # FIXME: get datadir from args
    data_dir = os.path.abspath(os.path.join(os.path.basename(__file__), 'HeatFlux'))

    heat_flux = np.loadtxt(os.path.join(data_dir, 'Antarctic_GHF.xyz'))
    heat_flux_unct = np.loadtxt(os.path.join(data_dir, 'Antarctic_GHF_uncertainty.xyz'))

    heat_flux_points = zip(heat_flux[:, 1], heat_flux[:, 0])
    new_bheatflx = scipy.interpolate.griddata(heat_flux_points, heat_flux[:, 2], (base.y_grid, base.x_grid), method='nearest')

    base.bheatflx = nc_base.creatVariable('bheatflx', 'f4', ('y','x'))
    base.bheatflx[:] = -new_bheatflx[:]  # invert sign for model convention
    base.bheatflx.grid_mapping = 'epsg_3031'
    base.bheatflx.coordinates = 'lon lat'

