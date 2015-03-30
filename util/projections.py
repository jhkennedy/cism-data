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


def greenland(args, lc_bamber, base):
    """The projections and tranformation grids for Greenland.

    This function creates the proj projections and a transformed DataGrid() for
    Greenland.

    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments.
    lc_bamber :
        Location of the Bamber dataset.
    base :
        A DataGrid() class instance that holds the base data grid information.

    Returns
    -------
    trans :
        A DataGrid() class instance that holds the base data grid transformed 
        to EPSG:3413.
    proj_eigen_gl04c :
        A proj class instance that holds the Bamber DEM projection.
    proj_epsg3413 :
        A proj class instance that holds the EPSG:3413 projection.
 
    """
    proj_epsg3413 = pyproj.Proj('+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84 +units=m')
        # InSAR data in this projections

    # EIGEN-GL04C referenced data:
    #----------------------------
    # unfortunately, bed, surface, and thickness data is referenced to 
    # EIGEN-GL04C which doesn't exist in proj4. However, EGM2008 should
    # be within ~1m everywhere (and within 10-20 cm in most places) so 
    # we use the egm08 projection which is available in proj4
    path_bamber = os.path.dirname(lc_bamber)
    if not ( os.path.exists(path_bamber+'/egm08_25.gtx') ):
        raise Exception("No "+path_bamber+"/egm08_25.gtx ! Get it here: http://download.osgeo.org/proj/vdatum/egm08_25/egm08_25.gtx") 
    
    #NOTE: Bamber projection appears to not actually have any fasle northings or eastings. 
    #proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +x_0=800000.0 +y_0=3400000.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')
    proj_eigen_gl04c = pyproj.Proj('+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0 +geoidgrids='+path_bamber+'/egm08_25.gtx')

    # transform meshes. 
    speak.verbose(args,"   Creating the transform meshes: base Bamber grid to EPSG-3413.")

    trans = DataGrid()
    trans.ny = base.ny
    trans.nx = base.nx
    trans.dims = base.dims

    trans.x_grid, trans.y_grid = pyproj.transform(proj_eigen_gl04c, proj_epsg3413, base.x_grid.flatten(), base.y_grid.flatten())
    trans.y_grid = trans.y_grid.reshape((base.ny,base.nx))
    trans.x_grid = trans.x_grid.reshape((base.ny,base.nx))

    return (trans, proj_epsg3413, proj_eigen_gl04c)
