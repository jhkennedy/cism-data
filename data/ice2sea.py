from util import speak

def apply_mask(args, nc_mask, nc_base ):
    """Apply Zurich mask to thk and make usrf.
    """
    base_thk  = nc_base.variables['thk']
    thk_data = base_thk[:,:]

    base_topg = nc_base.variables['topg']
    topg_data = base_topg[:,:]

    mask = nc_mask.variables['IceSheetMask']

    speak.verbose(args,"   Applying mask to thk.")
    thk_data = thk_data * mask[:,:]
    base_thk[:,:] = thk_data

    speak.verbose(args,"   Creating usrf.")
    base_usrf = nc_base.createVariable('usrf', 'f4',('y','x',))
    base_usrf[:,:] = thk_data + topg_data

