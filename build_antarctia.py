#!/usr/bin/env python2

"""
Build a CISM Dataset for Antarctica in the EPSG:3031 projection
"""

import os
import datetime
import subprocess
import argparse

from util import speak
from util import finalize
from util import projections
from util.ncfunc import get_nc_file
from util import custom_argparse_types as cats


