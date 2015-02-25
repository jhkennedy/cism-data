Scripts to builds CISM base data
================================

This is a collection of (mostly python) scripts to build a base dataset for the
Community Ice Sheet Model, [CISM](https://github.com/CISM). The requisite data
is not included. 

Instructions
------------
First, you will need to link or copy the data into the `data` directory or edit
the file paths on lines 22 through 30. Currently they are set to be: 

```python
#==== Data Locations ====
# Link data here or edit 
#========================
lc_bamber   = 'data/BamberDEM'           # file name: Greenland_bedrock_topography_V3.nc
lc_seaRise  = 'data/SeaRise'             # file name: Greenland1km.nc
lc_racmo2p0 = 'data/RACMO2.0'            # file name: Racmo2MeanSMB_1961-1990.nc
lc_InSAR    = 'data/InSAR/Joughin2012'   # file name: greenland_vel_mosaic500.nc #NOTE:  will build this file from mosaicOffsets.* files
lc_massCon  = 'data/IceBridge/Greenland' # file name: MCdataset-2014-11-19.nc
lc_mask     = 'data/Ice2Sea'             # file name: ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc
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

Some scripts were reused (with minor edits) and are in the `util/` directory while the rest live in `util/old/`. 
