"""
data.bedmap : Antarctic bed, surface, and thickness dataset.

This module provides ...

Functions list:
    * ...

Notes
-----

The data uses the EIGEN-GL04C projection:
    * Polar stereographic
    * WGS84 ellipsoid
    * EIGEN-GL04C geoidgrids
    * Standard parallel = -71 degrees
    * Latitude of projection origin = -90 degrees
    * Central meridian = 0 degrees
    * false eastings = 0
    * flase northings = 0
    * 1000 m postings with
        + lower-left corner y,x: -3333500.0,3333500.0 (m)
        + upper-right corner y,x: -3333500.0,3333500.0 (m)

The data is provided on a 6667 x 6667 grid.

References
----------
This dataset was provided by the authors of this publication:

Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, 2013: Ice-Shelf Melting
Around Antarctica. Science, 341, 266-270, doi:10.1126/science.1235798.

for optimization of validation purposes only. The authors must be contacted
prior to use in publications (as per authors request) and they should be
acknowledged as the source of the data in presentations.
"""

