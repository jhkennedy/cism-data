"""
uitl.projections : Utilities to ease pojection transformations.

This module provides classes and functions to help transform projections. 

Class list:
    * DataGrid()

Functions list:
    * greenland(args, lc_bamber, base)
"""

import os
import scipy
import pyproj

from util import speak

class DataGrid():
    """A class to hold data grids.
    """
    #FIXME: need a creation step here    
    def make_grid(self):
        """A way to make a basic grid from x,y data.
        """
        self.y_grid, self.x_grid = scipy.meshgrid(self.y[:], self.x[:], indexing='ij')
        self.dims = (self.ny, self.nx)
        self.dy = self.y[1]-self.y[0]
        self.dx = self.x[1]-self.x[0]

    def make_grid_flip_y(self):
        """A way to make a basic grid from x,y data, inverting y.
        """
        self.y_grid, self.x_grid = scipy.meshgrid(self.y[::-1], self.x[:], indexing='ij')
        self.dims = (self.ny, self.nx)
        self.dy = self.y[0]-self.y[1]
        self.dx = self.x[1]-self.x[0]


def grid_center_latlons(nc_base, base, proj):
    base.lon_grid = nc_base.createVariable('lon', 'f4', ('y','x',))
    base.lon_grid.long_name = 'grid center longitude'
    base.lon_grid.standard_name = 'longitude'
    base.lon_grid.units = 'degrees_east'
    base.lon_grid.grid_mapping = 'epsg_3413'
    base.lon_grid.note = 'Created by Joseph H. Kennedy using pyproj'

    base.lat_grid = nc_base.createVariable('lat', 'f4', ('y','x',))
    base.lat_grid.long_name = 'grid center latitude'
    base.lat_grid.standard_name = 'latitude'
    base.lat_grid.units = 'degrees_north'
    base.lat_grid.grid_mapping = 'epsg_3413'
    base.lat_grid.note = 'Created by Joseph H. Kennedy using pyproj'

    lon_grid, lat_grid = proj(base.x_grid.ravel(), base.y_grid.ravel(), inverse=True)
    lon_grid.shape = base.x_grid.shape 
    lat_grid.shape = base.x_grid.shape

    base.lon_grid[:,:] = lon_grid
    base.lat_grid[:,:] = lat_grid


def transform(base, proj1, proj2):
    trans = DataGrid()
    trans.ny = base.ny
    trans.nx = base.nx
    trans.dims = base.dims

    trans.x_grid, trans.y_grid = pyproj.transform(proj1, proj2, base.x_grid.flatten(), base.y_grid.flatten())
    trans.y_grid = trans.y_grid.reshape((base.ny,base.nx))
    trans.x_grid = trans.x_grid.reshape((base.ny,base.nx))

    return trans

def greenland():
    """The projections and tranformation grids for Greenland.

    This function creates the proj projections and a transformed DataGrid() for
    Greenland.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments.
    lc_bamber :
        Location of the Bamber dataset.

    Returns
    -------
    proj_eigen_gl04c :
        A proj class instance that holds the Bamber DEM projection.
    proj_epsg3413 :
        A proj class instance that holds the EPSG:3413 projection.
 
    """
    #NOTE: NSIDC sea ice polar stereographic north
    proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    #NOTE: WGS84 Arctic polar stereographic
    proj_epsg3995 = pyproj.Proj('+proj=stere +lat_0=90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

    # EIGEN-GL04C referenced data:
    #----------------------------
    # unfortunately, bed, surface, and thickness data is referenced to 
    # EIGEN-GL04C which doesn't exist in proj4. However, EGM2008 should
    # be within ~1m everywhere (and within 10-20 cm in most places) so 
    # we use the egm08 projection which is available in proj4
    path_egm08 = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'egm08_25.gtx')
    if not ( os.path.exists(path_egm08) ):
        raise Exception("No "+path_egm08+"! Get it here: http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx") 
    
    #NOTE: Bamber projection appears to not actually have any fasle northings or eastings. 
    #proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')
    proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +geoidgrids='+path_egm08)

    return (proj_epsg3413, proj_eigen_gl04c)

def equal_area(min_lat, max_lat, lon_0):
    proj_aea = pyproj.Proj('+proj=aea +lat_1='+str(min_lat)+' +lat_2='+str(max_lat)+' +lat_0='+str((max_lat+min_lat)/2.)+' +lon_0='+str(lon_0))
    return proj_aea


def antarctica():
    #NOTE: NSIDC sea ice polar stereographic south
    proj_epsg3412 = pyproj.Proj('+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
    #NOTE: WGS84 Antarctic polar stereographic
    proj_epsg3031 = pyproj.Proj('+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')

    # EIGEN-GL04C referenced data:
    #----------------------------
    # unfortunately, bed, surface, and thickness data is referenced to 
    # EIGEN-GL04C which doesn't exist in proj4. However, EGM2008 should
    # be within ~1m everywhere (and within 10-20 cm in most places) so 
    # we use the egm08 projection which is available in proj4
    path_egm08 = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'egm08_25.gtx')
    if not ( os.path.exists(path_egm08) ):
        raise Exception("No "+path_egm08+"! Get it here: http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx") 
    
    proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_0=-90 +lat_ts=71.0 +lon_0=0.0 +k_0=1.0 +geoidgrids='+path_egm08)
    
    return (proj_epsg3412, proj_eigen_gl04c)
