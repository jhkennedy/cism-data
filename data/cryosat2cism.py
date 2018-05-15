#!/usr/bin/env python

import numpy
import scipy.interpolate
import argparse

from netCDF4 import Dataset

from util import projections
from util import speak


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    volume = parser.add_mutually_exclusive_group()
    volume.add_argument("-v", "--verbose", help="Increase the output verbosity", action="store_true")
    volume.add_argument("-q", "--quiet",   help="Run silently",                  action="store_true")

    return parser.parse_args(args)


def main(args):
    # FIXME: from args
    lc_base    = 'Antarctica/antarctica_1km_2017_05_03.nc'
    lc_cryosat = 'cryosat2-dzdt/CS2_dzdt.nc'

    proj_epsg3412, proj_epsg3031 = projections.antarctica()

    # NOTE: Grid coordinates here relate to the Lower-Left corner of the grid cell.
    nc_cryosat = Dataset(lc_cryosat, 'r')

    cryosat = projections.DataGrid()
    cryosat.yll = nc_cryosat.variables['y1']
    cryosat.xll = nc_cryosat.variables['x1']
    cryosat.dy = cryosat.yll[1] - cryosat.yll[0]
    cryosat.dx = cryosat.xll[1] - cryosat.xll[0]
    cryosat.y = cryosat.yll[:] + cryosat.dy/2.
    cryosat.x = cryosat.xll[:] + cryosat.dx/2.
    cryosat.ny = cryosat.y[:].shape[0]
    cryosat.nx = cryosat.x[:].shape[0]
    cryosat.make_grid()

    cryosat.dzdt = numpy.ma.masked_values(nc_cryosat.variables['dzdt'][-1,:,:], -9999.)

    nc_base = Dataset(lc_base, 'r+')

    base = projections.DataGrid()
    base.y = nc_base.variables['y1']
    base.x = nc_base.variables['x1']
    base.dy = base.y[1]-base.y[0]
    base.dx = base.x[1]-base.x[0]
    base.ny = base.y[:].shape[0]
    base.nx = base.x[:].shape[0]
    base.N = base.ny*base.nx
    base.make_grid()

    cism2cry = projections.transform(base, proj_epsg3412, proj_epsg3031)
    cryosat_interp = scipy.interpolate.RectBivariateSpline(cryosat.y[:], cryosat.x[:], cryosat.dzdt, kx=1, ky=1, s=0)  # regular 2d linear interp. but faster

    DZDT = numpy.zeros( base.dims )
    for ii in range(0, base.nx):
        speak.progress(args, ii, base.nx)
        DZDT[:,ii] = cryosat_interp.ev(cism2cry.y_grid[:,ii], cism2cry.x_grid[:,ii])
    speak.progress(args, base.nx, base.nx)

    mask_y = (cism2cry.y_grid > cryosat.y[4:-4].max()) | (cism2cry.y_grid < cryosat.y[4:-4].min())
    mask_x = (cism2cry.x_grid > cryosat.x[4:-4].max()) | (cism2cry.x_grid < cryosat.x[4:-4].min())

    DZDT[mask_y | mask_x] = -9999
    
    
    base.dzdt = nc_base.createVariable('dzdt', 'f4', ('time','y1','x1',) )
    base.dzdt[0,:,:] = DZDT[:,:]
    base.dzdt.long_name = 'elevation rates of change'
    base.dzdt.units = 'meters/year'
    base.dzdt.missing_value = -9999
    base.dzdt.source = 'M. McMillan and A. Shepherd'
    base.dzdt.reference = 'Mcmillan, M., A. Shepherd, A. Sundal, K. Briggs, A. Muir, A. Ridout, A. Hogg, and D. Wingham, 2014: Increased ice losses from Antarctica detected by CryoSat-2. Geophys. Res. Lett, doi:10.1002/(ISSN)1944-8007.'
    base.dzdt.note = 'As per the request of the authors (M. McMillan & A. Shepherd), these data are not to be shared outside of this project (PISCEES). They are to be used for optimization and model validation purposes only, as the original authors still have plans to use them for other studies of their own. They are ok with us using them for optimization and validation with the caveat that we should communicate further with them about their use prior to publishing any stuides that use them. Also, if the data are used for presentations, we should acknowldege the authors as the source of the data. For any further questions, please check with S. Price or D. Martin.'
 
    nc_cryosat.close()
    nc_base.close()


if __name__ == "__main__":
    main(parse_args())
