Scripts to builds CISM base data
================================

This is a collection of (mostly python) scripts to build a base dataset for the
Community Ice Sheet Model, [CISM](https://github.com/CISM). The requisite data
is not included. 

Instructions
------------
First, you will need to link or copy the data into the `data` directory or edit
the file paths on lines 17 through 25. Currently they are set to be: 

```python
#==== Data Locations ====
# Link data here or edit 
#========================
lc_bamber   = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_mask     = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'
lc_massCon  = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'
lc_InSAR    = 'data/InSAR/Joughin2015/greenland_vel_mosaic500.nc' #NOTE:  will build this file from mosaicOffsets.* files
lc_racmo2p0 = 'data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc'
lc_seaRise  = 'data/SeaRise/Greenland1km.nc'
```

Once the data location paths are correct, you can build the datasets by running:
```bash
$ python build_greenland_datasets.py
```

which will take aproximitely 20 minutes to run. It will create a template netCDF4 
file called `templates/greenland_1km.mcb.nc` that has no time dimension. From
this file, a set of time-stamped `.mcb.nc` and `.mcb.config` files in the
`complete/` directory are created with an added time dimension for 1, 2, 4, and
8 km spacing over Greenland.

`build_greenland_datasets.py` also has some optional arguments:

```
usage: build_greenland_datasets.py [-h] [-v | -q]

optional arguments:
  -h, --help     show this help message and exit
  -v, --verbose  Increase the output verbosity
  -q, --quiet    Run silently
```

Data:
-----
Each `data/*.py` file's doc string contains a detailed description of the data. 

Authors:
--------
[Joseph H. Kennedy](https://github.com/jhkennedy) reworked scripts by 
[Steven F. Price](https://github.com/stephenprice), 
[Matt Hoffman](https://github.com/matthewhoffman), and 
[Matt Norman](https://github.com/matthewhoffman). 

The old scripts this package was built off of live in `util/old`.
