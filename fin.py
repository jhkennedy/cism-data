#!/usr/bin/env python2
# encoding: utf-8

"""
Build the antarctica dataset quick and dirty like based on an old version
"""

import shutil
import scipy.interpolate

import numpy as np

from argparse import Namespace
from datetime import datetime

from util import ncfunc
from util import projections
from util import finalize


f_1km = 'antarctica_1km_2018_04_13.nc'
f_template = ''
coarse_list = [2, 4, 5, 8]
args = Namespace()
args.quite = True
args.verbose = False


finalize.coarsen(args, 'epsg_3413', f_1km, f_template, coarse_list)

