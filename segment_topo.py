#!/usr/bin/env python2

from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndimage
from shapely.geometry import Point, shape 

from util import speak
from util import projections
from util.ncfunc import get_nc_file

RHO_OCEAN = 1027. # kg/m^3
RHO_ICE = 917.    # kg/m^3
R_RHO = RHO_ICE/RHO_OCEAN

def abs_existing_file(file):
    file = os.path.abspath(file)
    if not os.path.isfile(file):
        print("Error! File does not exist: \n    "+file)
        sys.exit(1)
    return file


def parse_args(args=None):
    parser = argparse.ArgumentParser(description=__doc__, 
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)   

    parser.add_argument('-i', '--input', type=abs_existing_file,
            default='templates/greenland_1km.epsg3413.nc',
            help='NetCDF dataset with topg variable to segment.')
    parser.add_argument('-a', '--altitude', type=float, default=-200.,
            help='Altitude (m. a.s.l.) at which to segment.')
    parser.add_argument('-w', '--write', action='store_true', 
            help='Write the mask to the input netCDF dataset.')
    parser.add_argument('-s', '--show', action='store_true', 
            help='Show the generated seglmentation plot.')

    volume = parser.add_mutually_exclusive_group()
    volume.add_argument('-v', '--verbose', action='store_true', 
            help='Increase the output verbosity')
    volume.add_argument('-q', '--quiet', action='store_true', 
            help='Run silently')

    return parser.parse_args(args)


def main(args):
    if args.write:
        nc_base = get_nc_file(args.input,'r+')
    else:
        nc_base = get_nc_file(args.input,'r')

    base = projections.DataGrid()
    base.y = nc_base.variables['y']
    base.x = nc_base.variables['x']
    base.ny = base.y[:].shape[0]
    base.nx = base.x[:].shape[0] 
    base.make_grid()


    base.topg = nc_base.variables['topg']
    base.thk = nc_base.variables['thk']

    topg = base.topg[::-1,:].filled()
    ice = base.thk[::-1,:]
    
    islands = np.greater(topg, args.altitude)
    segments, ids = ndimage.label(islands)
    #sgmt = ndimage.measurements.find_objects(segments)

    speak.verbose(args, '    Number of distinct shallow regions found: {}'.format(ids))

    mask = np.array(segments * -1.)
    # Set ocean
    mask[np.equal(mask,0)] = 1
    
    labeled = np.copy(mask) # For viewing

    #NOTE: nx = 2480, ny = 2975; -1 for index
    gl_x = [1200, 1990, 1990, 1625, 1320,  720,  440,  440, 359, 540, 586, 618, 710, 774, 819, 878, 969, 1200]
    gl_y = [  60,  344, 1680, 2225, 2940, 2940, 1995, 1168, 767, 523, 507, 493, 416, 391, 346, 305, 273,   60]
    shape_gl = shape({'type':'polygon', 'coordinates':[zip(gl_y, gl_x)]})

    ei_x = [960, 773, 140, 140, 158, 122, 144, 115,  84,  65,  45,  34,  44,  60,  76, 183, 207, 359, 540, 586, 618, 710, 774, 819, 878, 969, 960]
    ei_y = [ 84,  24,  24, 200, 290, 377, 431, 503, 505, 495, 501, 536, 558, 594, 652, 737, 855, 767, 523, 507, 493, 416, 391, 346, 305, 273,  84]
    shape_ei = shape({'type':'polygon', 'coordinates':[zip(ei_y, ei_x)]})

    
    for ii in range(mask.shape[0]):
        speak.progress(args, ii, mask.shape[0], width=60, char='=', indent=4)
        for jj in range(mask.shape[1]):
            if topg[ii,jj] == -9999.:
                mask[ii,jj] = 0
            
            elif Point((ii,jj)).within(shape_gl):
                if ice[ii,jj] > 0:
                    if ice[ii,jj] >= -1*topg[ii,jj]/R_RHO:
                        mask[ii,jj] = 3
                    else:
                        mask[ii,jj] = 4
                
                elif topg[ii,jj] > 0:
                    mask[ii,jj] = 2
                else:
                    mask[ii,jj] = 1
            
            elif mask[ii,jj] == 1:
                continue
            
            elif Point((ii,jj)).within(shape_ei):
                if topg[ii,jj] > 0:
                    mask[ii,jj] = -2
                else:
                    mask[ii,jj] = -1
            
            else:
                mask[ii,jj] = -3
    
    speak.progress(args, mask.shape[0], mask.shape[0], width=60, char='=', indent=4)

    if args.write:
        base.mask = nc_base.createVariable('mask', 'i4', ('y','x',) )
        base.mask[:,:] = mask[::-1,:]
        base.mask.long_name = 'location type mask'
        base.mask.grid_mapping = 'epsg_3413'
        base.mask.coordinates = 'lon lat'
        base.source = 'Joseph H. Kennedy, ORNL'
        base.comment = 'Mask values: -3, shallow ocean or land outside paleo domain; '+ \
                                    '-2, bare paleo land; '+ \
                                    '-1, shallow paleo ocean; '+ \
                                     '0, missing topg data; '+ \
                                     '1, ocean; '+ \
                                     '2, bare land; '+ \
                                     '3, grounded ice; '+ \
                                     '4, floating ice.'
                                    



    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize=(12,10), dpi=150)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
   
    ax1.imshow(topg, cmap='spectral', interpolation='nearest')
    ax1.set_title('Topography')

    ax2.imshow(islands, cmap='spectral', interpolation='nearest')
    ax2.set_title('Shallow regions')
    
    ax3.imshow(labeled, cmap='spectral', interpolation='nearest')
    ax3.set_title('Labeled regions')
    
    ax4.imshow(ice, cmap='spectral', interpolation='nearest')
    ax4.set_title('Ice thickness')
    
    ax5.imshow(mask, cmap='spectral', interpolation='nearest')
    ax5.set_title('Mask')

    ax6.imshow(mask, cmap='spectral', interpolation='nearest')
    ax6.autoscale(False)
    ax6.plot(gl_x, gl_y, '-ro')
    ax6.plot(ei_x, ei_y, '-bx')
    ax6.set_title('Mask with shapes')
   
    plt.tight_layout()
    plt.savefig('segments.png', bbox_inches='tight')
    
    if args.show:
        plt.show()

    nc_base.close()


if __name__ == '__main__':
    main(parse_args())

