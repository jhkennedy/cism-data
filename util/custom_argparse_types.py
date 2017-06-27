"""
A collection of functions which can be used as a custom type for argparse.
"""

import os
import sys
import errno

from argparse import ArgumentTypeError

def mkdir_p(path):
    """
    Make a directory, and parent directories as needed, with no error if
    directory already exists. This works like the Unix command `mkdir` with the
    `-p` option specified.
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


def abs_creation_dir(path):
    """
    Ensure the directory exists and return the absolute path to the directory.
    """
    path = os.path.abspath(path)
    mkdir_p(path)
    return path


def abs_existing_file(file_):
    """
    Ensure the file exists, it is a file, and return the absolute path to the
    directory. Else, raise an exception.
    """
    file_ = os.path.abspath(file_)
    if not os.path.isfile(file_):
        raise IOError(errno.ENOENT, 'File does not exist', file_)
    return file_
