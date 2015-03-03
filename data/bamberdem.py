import scipy

from util.ncfunc import copy_atts

def build_base (nc_bamber, nc_base, base):
    """Build the Bamber base grid.
    """
    bamber_y = nc_bamber.variables['projection_y_coordinate']
    bamber_ny = bamber_y[:].shape[0] # number of y points for 1km grid

    bamber_x = nc_bamber.variables['projection_x_coordinate']
    bamber_nx = bamber_x[:].shape[0] # number of x points for 1km grid

    # make bamber 1km grid for base
    base.ny = bamber_ny
    base.nx = bamber_nx

    # create our base dimensions
    nc_base.createDimension('y',base.ny)
    nc_base.createDimension('x',base.nx)

    # create some base variables
    base.y = nc_base.createVariable('y', 'f4', 'y')
    base.y[:] = bamber_y[:]
    copy_atts(bamber_y, base.y) #FIXME: units say km, but it's actuall in m

    base.x = nc_base.createVariable('x', 'f4', 'x')
    base.x[:] = bamber_x[:]
    copy_atts(bamber_x, base.x) #FIXME: units say km, but it's actuall in m

    # create some grids for interpolation
    base.y_grid, base.x_grid = scipy.meshgrid(base.y[:], base.x[:], indexing='ij')

    return base

