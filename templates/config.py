import os

def replace_dict(r_ewn, r_nsn, r_dew, r_dns, r_in, r_out, r_km):
    """Create replacements dictionary.
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
    """
    
    in_split = os.path.splitext(f_in)
    f_out = in_split[0]+'.out'+in_split[1]

    config_dict = replace_dict(base.ny, base.nx, km*1000, km*1000, f_in, f_out, str(km)+' km' )

    base_config = open('templates/'+f_template,'r')
    
    t_split = os.path.splitext( os.path.splitext(f_template)[0] )

    f_config = 'complete/'+t_split[0]+'_'+str(km)+'km'+t_split[1]+'.config'
    #NOTE: should make something like: complete/greenland_1km.mcb.config

    out_config  = open(f_config,'w')
    for line in base_config :
        for src, target in config_dict.iteritems() :
            line = line.replace(src, target)
        out_config.write(line)

    base_config.close()
    out_config.close()
    os.chmod(f_config, 0o755)   # uses an octal number!


