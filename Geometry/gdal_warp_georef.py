from pathlib import Path
from osgeo import gdal, osr

# Adapted from https://svn.osgeo.org/gdal/trunk/autotest/alg/warp.py
def warp_with_gcps(input_path, output_path, gcps, gcp_epsg=4326, output_epsg=4326):
    # Open the source dataset and add GCPs to it
    src_ds = input_path#gdal.OpenShared(str(input_path), gdal.GA_Update)
    gcp_srs = osr.SpatialReference()
    gcp_srs.ImportFromEPSG(gcp_epsg)
    gcp_crs_wkt = gcp_srs.ExportToWkt()
    src_ds.SetGCPs(gcps, gcp_crs_wkt)

    # Define target SRS
    dst_srs = osr.SpatialReference()
    dst_srs.ImportFromEPSG(output_epsg)
    dst_wkt = dst_srs.ExportToWkt()

    error_threshold = 0.125  # error threshold --> use same value as in gdalwarp
    resampling = gdal.GRA_Bilinear

    # Call AutoCreateWarpedVRT() to fetch default values for target raster dimensions and geotransform
    tmp_ds = gdal.AutoCreateWarpedVRT(src_ds,
                                      None,  # src_wkt : left to default value --> will use the one from source
                                      dst_wkt,
                                      resampling,
                                      error_threshold,
                                      )
    dst_xsize = tmp_ds.RasterXSize
    dst_ysize = tmp_ds.RasterYSize
    dst_gt = tmp_ds.GetGeoTransform()
    tmp_ds = None

    # Now create the true target dataset
    dst_path = str(Path(output_path).with_suffix(".tif"))
    dst_ds = gdal.GetDriverByName('GTiff').Create(dst_path, dst_xsize, dst_ysize, src_ds.RasterCount)
    dst_ds.SetProjection(dst_wkt)
    dst_ds.SetGeoTransform(dst_gt)
    dst_ds.GetRasterBand(1).SetNoDataValue(0)

    # And run the reprojection
    gdal.ReprojectImage(src_ds,
                        dst_ds,
                        None,  # src_wkt : left to default value --> will use the one from source
                        None,  # dst_wkt : left to default value --> will use the one from destination
                        resampling,
                        0,  # WarpMemoryLimit : left to default value
                        error_threshold,
                        None,  # Progress callback : could be left to None or unspecified for silent progress
                        None)  # Progress callback user data
    dst_ds = None
# input_path = Path("x.tif")
# output_path = Path("y.tif")
# # GCP input
# xyz = [...]
# row_col = [...]
#
# gcps = []
# for (x, y, z), (row, col) in zip(xyz, row_col):
#     gcps.append(gdal.GCP(x, y, z, col, row))
#
# warp_with_gcps(input_path, output_path, gcps, gcp_epsg=3301, output_epsg=3301)
