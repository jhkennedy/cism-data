import numpy as np
import scipy

def check_missing(data, anchor, i_s, j_s):
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
    
    w = np.zeros(4)

    w[0] = (1.-alpha)*(1.-beta)
    w[1] = (1.-alpha)*(   beta)
    w[2] = (   alpha)*(   beta)
    w[3] = (   alpha)*(1.-beta)

    return w
