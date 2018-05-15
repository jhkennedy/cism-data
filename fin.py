#!/usr/bin/env python2
# encoding: utf-8

"""
Build the antarctica dataset quick and dirty like based on an old version
"""


import argparse

from util import finalize

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--coarse-list', type=int, nargs='+',
        help='Coarsen the 1 km dataset to this resolution (in km).')
args = parser.parse_args()
args.quite = True
args.verbose = False

f_1km = 'antarctica_1km_2018_05_14.nc'
f_template = ''

finalize.coarsen(args, 'epsg_3031', f_1km, f_template, args.coarse_list)

