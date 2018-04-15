"""
data.subself : Antarctic sub-shelf marine melt rates import module.

This module provides ...

Functions list:
    * ...

Notes
-----

The data uses the EPSG:3031 projection:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = -71 degrees
    * Latitude of projection origin = -90 degrees
    * Central meridian = 0 degrees
    * false eastings = 0
    * flase northings = 0
    * 1000 m postings with
        + lower-left corner y,x: -2800000.0,-2800000.0 (m)
        + upper-right corner y,x: -2800000.0,-2800000.0 (m)

The data is provided on a 5601 x 5601 grid and includes latitude and longitude
values.

Importantly, the y axis in this dataset is opposite the standard CISM convention
and will need to be flipped.

References
----------
This dataset was provided by the authors of this publication:

Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, 2013: Ice-Shelf Melting
Around Antarctica. Science, 341, 266-270, doi:10.1126/science.1235798.

for optimization of validation purposes only. The authors must be contacted
prior to use in publications (as per authors request) and they should be
acknowledged as the source of the data in presentations.
"""

