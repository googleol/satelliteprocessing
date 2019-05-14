# satelliteprocessing

Sorry that the code is very dirty.
This was a R&D excersise just to test the algorithmic workflow of generation of
a satellite image with its orbit and attitude data.

The main terrain correction ( intersection of look vector with WGS84+DEM )is in the FUNCTION
RAY_INTERSECT_WGS84_TERRAIN .

The general work flow to generate the Look ray from payload angles is all in the CCD_TO_GROUND_TERRAIN file .

Important parameters here are satellite location in orbit and it's roll pitch yaw values.

All the values in the config file are fake. The input is a h5 file on my side. You have to replace
how the input(data and orbit attitude parameters ) are given to the program in the CCD_TO_GROUND_TERRAIN file.

Import functions which are useful are in rad mat utils.

gdal_warp_georef is used for final generation of a geo tiff from the generated grid.
