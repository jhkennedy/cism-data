import os

def replace_dict(r_nsn, r_ewn, r_dew, r_dns, r_in, r_out, r_km):
    """Create replacements dictionary.

    This function creates a replacement dictionary for replacing key words in 
    the config file template.

    Parameters
    ----------
    r_ewn :
       Number of grid points in the y direction.
    r_nsn :
       Number of grid points in the x direction.
    r_dew :
        Grid separation in the y direction (in meters). 
    r_dns :
        Grid separation in the x direction (in meters). 
    f_in :
       The input dataset name.
    f_out :
        The output dataset name (created by CISM).
    r_km :
        Grid resolution in kilometers.
    """
    config_dict= {'REPLACE_EWN' : str(r_ewn), 
                  'REPLACE_NSN' : str(r_nsn), 
                  'REPLACE_DEW' : str(r_dew), 
                  'REPLACE_DNS' : str(r_dns), 
                  'REPLACE_NAME': str(r_in ), 
                  'REPLACE_OUT' : str(r_out), 
                  'REPLACE_KM'  : str(r_km ) }
    return config_dict


def write(f_in, f_template, base, km):
    """write config file.

    This function creates a config file from the config file template.

    Parameters
    ----------
    f_in :
       The input dataset name.
    f_template :
        The name of the config file template file.
    base :
        A DataGrid() class holding the base data grid.
    km :
        Grid resolution in kilometers.
    """
    
    in_bpath, in_bfile = os.path.split(f_in)
    in_bname, in_bext = os.path.splitext(in_bfile)
    
    f_out = os.path.join(in_bpath, in_bname+'.out'+in_bext)

    config_dict = replace_dict(base.ny, base.nx, km*1000, km*1000, f_in, f_out, str(km)+' km' )

    base_config = open(f_template,'r')
    
    t_split = os.path.basename(f_template).split('.') 

    f_config = os.path.join(in_bpath, t_split[0]+'_'+str(km)+'km.'+'.'.join(t_split[1:]))
    #NOTE: should make something like: complete/greenland_1km.mcb.config

    out_config  = open(f_config,'w')
    for line in base_config :
        for src, target in config_dict.iteritems() :
            line = line.replace(src, target)
        out_config.write(line)

    base_config.close()
    out_config.close()
    os.chmod(f_config, 0o755)   # uses an octal number!


