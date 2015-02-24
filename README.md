Scripts to builds CISM base data
================================

This is a collection of (mostly python) scripts to build a base dataset for the
Community Ice Sheet Model, [CISM](https://github.com/CISM). The requisite data
is not included. 

Instructions
------------
First, you will need to link or copy the data into the `data` directory or edit
the file paths on lines 22 through 32. Currently they are set to be: 

```python
# load in datasets
nc_bamber   = Dataset( 'data/BamberDEM/Greenland_bedrock_topography_V3.nc', 'r') 
nc_seaRise  = Dataset( 'data/SeaRise/Greenland1km.nc', 'r')
nc_racmo2p0 = Dataset( 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc', 'r')

if not ( os.path.exists('data/InSAR/Joughin2012/greenland_vel_mosaic500.nc') ):
    subprocess.call("python util/convert_velocities.py", shell=True)
nc_insar    = Dataset( 'data/InSAR/Joughin2012/greenland_vel_mosaic500.nc' , 'r')

nc_massCon  = Dataset('data/IceBridge/Greenland/MCdataset-2014-11-19.nc','r')
nc_mask     = Dataset( 'data/Ice2Sea/greenland_geometry_icesheet_mask_Zurich.nc' , 'r'  )
```


`build_greenland_datasets.py` is the main build script. It is, for the most
part, one long script that takes a bit of time to run. You can run it from
within the top-level repository directory as such:

```bash
$ python build_greenland_datasets.py
```

It will create a netCDF4 file called `greenland_1km.mcb.nc` and a set of
time-stamped `.mcb.nc` and `.mcb.config` files in the `complete/` directory for
1, 2, 4, and 8 km spacing over Greenland. 

Authors:
--------
[Joseph H. Kennedy](https://github.com/jhkennedy), reworked scripts by 
[Steven F. Price](https://github.com/stephenprice), 
[Matt Hoffman](https://github.com/matthewhoffman), and 
[Matt Norman](https://github.com/matthewhoffman). 

Some scripts were reused are in the `util/` directory while the rest live in `util/old/`. 
