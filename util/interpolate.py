"""
uitl.interpolate : A set of interpolation helper functions.

This module provides functions to help create interpolation routines. 

Functions list:
    * check_missing(data, anchor, i_s, j_s)
    * linear_weights(alpha, beta)
"""

import numpy as np
import scipy

def check_missing(data, anchor, i_s, j_s):
    """Check for missing data.

    This function checks a masked data array around an anchor point for 
    missing values.

    Notes
    -----

    Points are ordered as such:
        0: (ii    , jj    )
        1: (ii    , jj+j_s)
        2: (ii+i_s, jj+j_s)
        3: (ii+i_s, jj    )

        Which, for an upper-right quadrent looks like:

        3 ---- 2
        |      |
        |      |
        0 ---- 1

        Numbering in other quadrents is the reflection through the axis or axes 
        with negitive skip values (i_s or j_s).

    Parameters
    ----------
    data :
        A masked numpy data array. 
    anchor :
        A tuple with the ii and jj index of the anchor point.
    i_s :
        A +/- 1 value which determines ii index of the point quadrent to
        inspect around the anchor point.
    j_s :
        A +/- 1 value which determines jj index of the point quadrent to
        inspect around the anchor point.

    Returns
    -------
    missing_points :
        A list of tuples containing the ii,jj indicies of points with missing
        data.
    interp_dict :
        A dictionary containg a list of tuples contining the ii,jj indicies of
        the checked quadrent's points.
    """
    missing_points = []
    interp_dict = {}

    dims_list = [(anchor[0],    anchor[1]    ),
                 (anchor[0],    anchor[1]+j_s),
                 (anchor[0]+i_s,anchor[1]+j_s),
                 (anchor[0]+i_s,anchor[1]    )]

    for dims in dims_list :
        if not any( d < 0 for d in dims) and data.mask[dims]:
            missing_points.append(dims)
        
        interp_dict.update({dims: data[dims]} )


    return (missing_points, interp_dict)


def linear_weights(alpha, beta):
    """Check for missing data.

    This function returns a list of linear interpolation weights.

    Parameters
    ----------
    alpha :
        Relitive, normalized location, in ii, of point to interpolate within a 
        points quadrent. 
    beta :
        Relitive, normalized location, in jj, of point to interpolate within a 
        points quadrent. 

    Returns
    -------
    w :
        A list of linear interpolation weights for an interpoation quadrent
        ordered like:
        0: (ii    , jj    )
        1: (ii    , jj+j_s)
        2: (ii+i_s, jj+j_s)
        3: (ii+i_s, jj    )

        Which, for an upper-right quadrent looks like:

        3 ---- 2
        |      |
        |      |
        0 ---- 1

        Numbering in other quadrents is the reflection through the axis or axes 
        with negitive skip values (i_s or j_s).
    """
    
    w = np.zeros(4)

    w[0] = (1.-alpha)*(1.-beta)
    w[1] = (1.-alpha)*(   beta)
    w[2] = (   alpha)*(   beta)
    w[3] = (   alpha)*(1.-beta)

    return w
