# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

DEM TILE PRODUCTS LIBRARY:
Library of Python functions to process various DEM tile products:
- SRTMGL1 V003 (SRTM 1 arc-sec) | NASA/LPDAAC
- ASTGTM2 V002 (ASTER GDEMv2 1 arc-sec) | NASA/LPDAAC
- ArcticDEM mosaic release 7 (2m, 10m, 32m) | NGA-NSF/PGC
- REMA mosaic release 1.1 (8m) | PGC
- AW3D30 (ALOS World 3D 30-m DEM) | JAXA
- SRTMv4.1 (SRTM 90m void-filled) | NASA/CGIAR-CSI
- TDX-90m (TanDEM-X 90m DEM) | DLR/ADS
- NASADEM (SRTM revisited) | NASA <-- still waiting for final product ; TODO: NASADEM to add

Processing steps description:
1/ Reconcile product naming convention with input
2/ Product extraction
3/ DEM "basic" filtering options using provided mask/err/num/metadata
4/ DEM warp (projection, resolution, extent, nodata_out, resampling method)
5/ Vertical referencing (geoid or ellipsoid)
For some data-sets, possible options (not by default):
6/ DEM sorting and merging to the extent of a 1x1° lat/lon tile instead of native tile type (ArcticDEM, REMA)
7/ DEM clipping to 1x1° lat/lon tile if native tile type is different (SRTMv4.1, ArcticDEM, REMA)

Generic input description:

- dir: path to data directory !!CHECK REQUIRED DIRECTORY ARCHITECTURE IN FUNCTION DESCRIPTION!!
usually either structure native to the bulk downloading (ArcticDEM, REMA, SRTMv4.1) or all archives at the root of pointed directory (SRTMGL1 & GL1N, ASTGTM2, AW3D30, TDX)

- tile_name: tile name (SRTMGL1 convention : South West corner of tile), e.g. N52W084

- dem_out: path to output file

format_out, tgt_EPSG, tgt_res, nodata_out, interp_method: warping parameters (GDAL naming)
in details with examples,
- format_out: GDAL format naming at: https://www.gdal.org/formats_list.html, e.g. 'GTiff', 'HDF4'
- tgt_EPSG: EPSG as int, e.g. 4326, 3413
- tgt_res= [xres,yres], same as input for "gdalwarp -tr", e.g. [30,-30] for 30m in UTM projection
- nodata_out= no-data value as int, e.g. -9999, -32768
- interp_method: same naming than for "gdalwarp -r", e.g. 'near', 'bilinear', 'cubic', 'cubicspline', etc...

- geoid: True or False, respectively set output vertical reference to either geoid or ellipsoid (input reference is detected)

- bit_params or filt_params: basic filtering parameters for num/mask file: see "rastlib" functions "filter_nanarray" and "mask_bitarray" for an example

and for some datasets,
- tag_latlon_tile: True or False, tile naming to interpret: respectively SRTMGL1 or native to data-set (e.g. ArcticDEM/REMA mosaic tile naming)

- tag_merge: True or False, merge to extent or not (if different than native)

- tag_clip: True or False, clip to extent or not (if different than native)
"""

import numpy as np
import os, shutil, sys
import zipfile
from misclib import SRTMGL1_naming_to_latlon
from shlib import create_tmp_dir_for_outfile, remove_tmp_dir_for_outfile, extract_file_from_tar_gz, extract_file_from_zip
from rastlib import read_nanarray, update_nanarray, filter_nanarray, mask_bitarray, warp_defaultUTM, clip_rast_to_extent, merge_rast_list
from vectlib import list_shp_field_inters_extent
from demlib import geoid_to_ellipsoid, ellipsoid_to_geoid

def AW3D30_v2_1_tile(dir_ALOS,tile_name,dem_out,bit_params=None,format_out='GTiff',tgt_EPSG=4326,tgt_res=None,nodata_out=-9999,interp_method=None,geoid=True):

    """
    :param dir_ALOS: path to directory directly containing ALOS World 3D 30-m 5x5° tar.gz archives
    :param tile_name: 1x1° tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param bit_params: bit filtering for AW3D30 MSK file ; see rastlib.mask_bitarray function for an example
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :return

    ALOS World 3D-30m (AW3D30) Version 2.1 from JAXA at: https://www.eorc.jaxa.jp/ALOS/en/aw3d30/aw3d30v21_format_e.pdf

    MSK file:
    0000: Valid
    0001: Cloud and show mask (invalid)
    0010: Land water and low correlation mask of the 5 m resolution DSM (valid)
    0011: Sea mask (valid)
    0100: National Land Numerical Information 10 m DEM (by Geographical Survey Institute of Japan) (valid)
    1000: Shuttle Radar Topography Mission (SRTM) SRTM-1 Version 3 （valid）
    1100: PRISM DSM (valid)

    STK, HDR, QAI, LST files not used in this function, could be extracted same manner as other in 2/

    OPTIMAL DIRECTORY ARCHITECTURE: all 5x5 .tar.gz archives in the same directory
    """

    #1/ LOCATE 1x1° TILE

    #ALOS AW3D30 v2.1 5x5 archive naming convention
    def latlon_to_ALOS_5x5tile_naming(lat,lon):

        lat_5=np.floor(lat/5)*5
        lon_5=np.floor(lon/5)*5

        if lat_5>=0:
            str_lat='N'
        else:
            str_lat='S'
        if lon_5>=0:
            str_lon='E'
        else:
            str_lon='W'

        tile_name_5 = str_lat + str(lat_5).zfill(3) + str_lon + str(lon_5).zfill(3)

        return tile_name_5

    #identify in which 5x5 tar.gz the 1x1 tile is located

    lat_tile,lon_tile=SRTMGL1_naming_to_latlon(tile_name)

    name_5x5 = latlon_to_ALOS_5x5tile_naming(lat_tile,lon_tile)

    list_5x5 = os.listdir(dir_ALOS)
    list_5x5_prefix = [l[0:8] for l in list_5x5]

    if name_5x5 in list_5x5_prefix:
        tar_5x5=list_5x5[list_5x5_prefix.index(name_5x5)]
    else:
        sys.exit('Could not find an ALOS AW3D30 5x5° archive containing the tile '+tile_name+ ' in the directory '+dir_ALOS)

    #2/ EXTRACT TILE
    tmp_dir=create_tmp_dir_for_outfile(dem_out)

    alos_tile_name=tile_name[0:1]+str(0)+tile_name[1:]  #SRTMGL1 naming convention to ALOS 1x1 naming convention (3 digits for latitude)

    #extract files of interest
    dsm_file=alos_tile_name+'_AVE_DSM.tif'
    msk_file=alos_tile_name+'_AVE_MSK.tif'

    tmp_dsm=os.path.join(tmp_dir,dsm_file)
    tmp_msk=os.path.join(tmp_dir,msk_file)

    extract_file_from_tar_gz(os.path.join(dir_ALOS,tar_5x5),dsm_file,tmp_dsm)
    extract_file_from_tar_gz(os.path.join(dir_ALOS,tar_5x5),msk_file,tmp_msk)

    #3/ FILTER TILE

    mask=read_nanarray(tmp_msk)

    #no filtering by default
    if bit_params is not None:

        mask_goodval=mask_bitarray(mask,bit_params[0],bit_params[1])
        dsm = read_nanarray(tmp_dsm)

        dsm_filtered= np.array(dsm)
        dsm_filtered[mask_goodval]=np.NaN

        update_nanarray(tmp_dsm,dsm_filtered)

    #4/ REPROJECT TILE:

    #raw data is GeoTiff, 4326, 1 arc-sec and -9999 nodata_out
    if format_out=='GTiff' and tgt_EPSG==4326 and tgt_res is None and nodata_out is -9999:
        tmp_dsm_proj = tmp_dsm
    else:
        tmp_dsm_proj = os.path.join(tmp_dir,alos_tile_name+'_proj.tif')
        warp_defaultUTM(tmp_dsm,tmp_dsm_proj,format_out,4326,tgt_EPSG,tgt_res,nodata_out,interp_method)

    #5/ ELLIPSOID OR GEOID

    #raw data is geoid EGM96
    if not geoid:
        geoid_to_ellipsoid(tmp_dsm_proj,dem_out)
    else:
        shutil.move(tmp_dsm_proj,dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def TDX_90m_tile(dir_TDX,tile_name,dem_out,filter_params=None,format_out='GTiff',tgt_EPSG=4326,tgt_res=None,nodata_out=-32767,interp_method=None,geoid=False):

    """
    :param dir_TDX: path to directory directly containing TanDEM-X zip archives
    :param tile_name: 1x1° tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param filter_params: filtering with TDX HEM file ; see rastlib.filter_nanarray function for an example
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :return:

    TanDEM-X 90m product from DLR at: https://geoservice.dlr.de/web/dataguide/tdm90/pdfs/TD-GS-PS-0021_DEM-Product-Specification.pdf

    Only HEM (Height Error Map) is extracted here
    Using the same approach with binary/nanarray filters, 2/ can extract and 3/ can process the AMP, AM2, WAM, COV, COM and LSM files

    OPTIMAL DIRECTORY ARCHITECTURE: all .zip archives in the same directory
    """

    # 1/ LOCATE 1x1° TILE
    def latlon_to_TDX_tile_naming(lat,lon):

        if lat>=80:
            lon_sw=np.floor(lon/4)*4
        elif lat >=60:
            lon_sw = np.floor(lon / 2) * 2
        else:
            lon_sw = lon

        lat_sw = lat

        if lat_sw>=0:
            str_lat='N'
        else:
            str_lat='S'
        if lon_sw>=0:
            str_lon='E'
        else:
            str_lon='W'

        tile_name_tdx= str_lat + str(int(abs(lat_sw))).zfill(2) + str_lon + str(int(abs(lon_sw))).zfill(3)

        return tile_name_tdx

    #identify in which TDX tile the 1x1 tile is located
    lat_tile,lon_tile=SRTMGL1_naming_to_latlon(tile_name)
    name_tdx = latlon_to_TDX_tile_naming(lat_tile,lon_tile)

    list_tdx = os.listdir(dir_TDX)
    list_prefix = [l[13:20] for l in list_tdx]

    if name_tdx in list_prefix:
        pass
        # zip_tdx=list_prefix[list_prefix.index(name_tdx)]
    else:
        sys.exit('Could not find a TDX tile archive containing the tile '+name_tdx+ ' in the directory '+dir_TDX)

    tile_name = name_tdx
    # nomenclature of TDX 3 arc-sec files
    tile_zip = os.path.join(dir_TDX,'TDM1_DEM__30_'+name_tdx +'.zip')

    dem_out = dem_out + '_' + name_tdx +'.tif'

    # 2/ EXTRACT TILE

    tmp_dir = create_tmp_dir_for_outfile(dem_out)

    dem_file='TDM1_DEM__30_'+tile_name+'_DEM.tif'
    hem_file='TDM1_DEM__30_'+tile_name+'_HEM.tif'

    tmp_dem = os.path.join(tmp_dir,dem_file)
    tmp_hem = os.path.join(tmp_dir,hem_file)

    extract_file_from_zip(tile_zip,dem_file,tmp_dem)
    extract_file_from_zip(tile_zip,hem_file,tmp_hem)

    # 3/ FILTER TILE
    if filter_params is not None:
        hem = read_nanarray(tmp_hem)
        _, filt = filter_nanarray(hem, filter_params[0], filter_params[1], filter_params[2])

        dem = read_nanarray(tmp_dem)
        dem_filtered = np.array(dem)
        dem_filtered[filt] = np.NaN

        update_nanarray(tmp_dem, dem_filtered)

    # 4/ REPROJECT TILE

    # raw data is GTiff, 4326, 3 arc-sec and -32768 nodata_out
    if format_out == 'GTiff' and tgt_EPSG == 4326 and tgt_res is None and nodata_out is -32767:
        tmp_dem_proj = tmp_dem
    else:
        tmp_dem_proj = os.path.join(tmp_dir, tile_name + '_proj.tif')
        warp_defaultUTM(tmp_dem, tmp_dem_proj, format_out=format_out, src_EPSG=4326, tgt_EPSG=tgt_EPSG, tgt_res=tgt_res, nodata_in=-32767, nodata_out=nodata_out, interp_method=interp_method)

    # 5/ ELLIPSOID OR GEOID

    # raw data is ellipsoid WGS84
    if geoid:
        ellipsoid_to_geoid(tmp_dem_proj, dem_out)
    else:
        shutil.move(tmp_dem_proj, dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def SRTMv4_1_tile(dir_SRTMv4_1,tile_name,dem_out,tag_clip=False,filter_params=None,format_out='GTiff',tgt_EPSG=4326,tgt_res=None,nodata_out=-32768,interp_method=None,geoid=True):

    """
    :param dir_SRTMv4_1: path to parent directory containing "6_5x5_TIFS" directory of zip archives (bulk downloading)
    :param tile_name: 1x1° tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param tag_clip: True to clip the 5x5° tile to the 1x1° extent of tile_name
    :param filter_params: no filtering for this product
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :return:

    SRTM 90m Digital Elevation Database v4.1 by NASA/CGIARCSI at: https://cgiarcsi.community/data/srtm-90m-digital-elevation-database-v4-1/

    Product Description temporarily not available (web issue)

    OPTIMAL DIRECTORY ARCHITECTURE: bulk download format, point to the directoy containing 6_5x5_TIFS directory
    """

    def latlon_to_SRTMv4_1_5x5tile_naming(lat,lon):

        ind_lon = np.floor(lon +180 / 5) +1
        ind_lat = np.floor(lat +90/ 5) +1

        return ind_lat, ind_lon

    # 1/ LOCATE 5x5° tile

    #read lat, lon according to SRTMv4.1 naming convention
    lat_tile, lon_tile = SRTMGL1_naming_to_latlon(tile_name)
    ind_lat_4_1, ind_lon_4_1 = latlon_to_SRTMv4_1_5x5tile_naming(lat_tile, lon_tile)

    # nomenclature of SRTMv4.1 archives
    tile_zip = os.path.join(dir_SRTMv4_1, '6_5x5_TIFS','srtm_' + str(ind_lon_4_1).zfill(2)+'_'+str(ind_lat_4_1).zfill(2) + '.zip')

    # 2/ EXTRACT TILE

    tmp_dir = create_tmp_dir_for_outfile(dem_out)

    dem_file = 'srtm_' + str(ind_lon_4_1).zfill(2)+'_'+str(ind_lat_4_1).zfill(2) + '.tif'

    tmp_dem = os.path.join(tmp_dir, dem_file)

    extract_file_from_zip(tile_zip, dem_file, tmp_dem)

    #2.5/ CLIP TO TILE EXTENT

    if not tag_clip:
        tmp_dem_clipped=os.path.join(tmp_dir,'srtm_' + str(ind_lon_4_1).zfill(2)+'_'+str(ind_lat_4_1).zfill(2)+'_clipped.tif')
        clip_rast_to_extent(tmp_dem,tmp_dem_clipped,[lat_tile,lon_tile,lat_tile+1,lon_tile+1],4326)
    else:
        tmp_dem_clipped=tmp_dem

    # 3/ FILTER TILE
    if filter_params is not None:
        sys.exit('No filter pre-defined for this DEM product.')

    # 4/ REPROJECT TILE

    # raw data is GeoTiff, 4326, 3 arc-sec and -32768 nodata_out
    if format_out == 'GTiff' and tgt_EPSG == 4326 and tgt_res is None and nodata_out is -32768 and tag_clip:
        tmp_dem_proj = tmp_dem_clipped
    else:
        tmp_dem_proj = os.path.join(tmp_dir, tile_name + '_proj.tif')
        warp_defaultUTM(tmp_dem_clipped, tmp_dem_proj, format_out, 4326, tgt_EPSG, tgt_res, nodata_out, interp_method)

    # 5/ ELLIPSOID OR GEOID

    # raw data is geoid EGM96
    if not geoid:
        geoid_to_ellipsoid(tmp_dem_proj, dem_out)
    else:
        shutil.move(tmp_dem_proj, dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def ArcticDEM_mosaic_r7_tile(dir_ArcDEM,tile_name,dem_out,filter_params=None,format_out='GTiff',tgt_EPSG=3413,tgt_res=None,nodata_out=-9999,interp_method=None,geoid=False,tag_lonlat_tile=False,path_tile_index=None,tag_merge=False,tag_clip=False):

    """
    :param dir_ArcDEM: path to parent directory "2m", "10m" or "32m" containing subdirectories of tar.gz archives (native FTP architecture)
    :param tile_name: either ArcticDEM tile name or 1x1° lat/lon tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param filter_params: no filtering for this product
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :param tag_lonlat_tile: True if tile_name follows SRTMGL1 tile naming, False if tile_name follows ArcticDEM tile naming
    :param path_tile_index: if tile_name is ArcticDEM format, specify path to ESRI ArcticDEM Tile Index
    :param tag_merge: if tile_name is ArcticDEM format, True to merge all ArcticDEM tiles to the 1x1° lat/lon extent
    :param tag_clip: if tile_name is ArcticDEM format, True to clip the 5x5° tile to the 1x1° lat/lon extent of tile_name
    :return:
    ArcticDEM release 7 product: ref:https://www.pgc.umn.edu/data/arcticdem/

    Processing should work for 2m, 10m and 30m versions of the mosaic
    (100m, 500m and 1km versions are bundled in one .tif file)

    Tile name and processing is ArcticDEM tile naming convention by default
    Provide path to ESRI tile index file to use 1x1° lat/lon tiles and SRTMGL1 naming convention

    No filtering is defined.

    OPTIMAL DIRECTORY ARCHITECTURE: point to "2m", "10m" or "32m" folder of similar architecture than: ftp://ftp.data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0
    """

    # 1/ LOCATE TILE

    if not tag_lonlat_tile:
        subtile_dir=os.path.join(dir_ArcDEM,tile_name)
        tile_tar_gz_list=[os.path.join(subtile_dir,tar_file) for tar_file in os.listdir(subtile_dir) if tar_file.endswith('.tar.gz')]
    else:
        lat_tile,lon_tile=SRTMGL1_naming_to_latlon(tile_name)
        extent=[lat_tile,lon_tile,lat_tile+1,lon_tile+1]
        #feature name in ArcticDEM_Tile_Index_Rel7
        feat_name='tile'
        subtile_name_list=list_shp_field_inters_extent(path_tile_index, feat_name, extent,4326)
        subtile_dir_list = [os.path.join(dir_ArcDEM,tile) for tile in subtile_name_list]
        tile_tar_gz_list=[]
        for i in range(len(subtile_dir_list)):
            tile_tar_gz_list=tile_tar_gz_list+[os.path.join(subtile_dir_list[i],tar_file) for tar_file in os.listdir(subtile_dir_list[i]) if tar_file.endswith('.tar.gz')]

    # 2/ EXTRACT TILE

    tmp_dir = create_tmp_dir_for_outfile(dem_out)

    list_tmp_dem = [os.path.join(tmp_dir, os.path.splitext(os.path.basename(tile_tar_gz))[0]+'_reg_dem.tif') for tile_tar_gz in tile_tar_gz_list]
    for tile_tar_gz in tile_tar_gz_list:
        extract_file_from_tar_gz(tile_tar_gz,os.path.splitext(os.path.basename(tile_tar_gz))[0]+'_reg_dem.tif',list_tmp_dem[tile_tar_gz_list.index(tile_tar_gz)])

    list_tmp_dem_tomerge=[]
    for tmp_dem in list_tmp_dem:
        # 3/ FILTER TILE
        if filter_params is not None:
            sys.exit('No filter pre-defined for this DEM product.')

        # 4/ REPROJECT TILE

        # raw data is GeoTiff, 3413, 1 arc-sec and -9999 nodata_out
        if format_out == 'GTiff' and tgt_EPSG == 3413 and tgt_res is None and nodata_out is -9999:
            tmp_dem_proj = tmp_dem
        else:
            tmp_dem_proj = os.path.join(tmp_dir, os.path.splitext(os.path.basename(tmp_dem))[0] + '_proj.tif')
            warp_defaultUTM(tmp_dem, tmp_dem_proj, format_out, 3413, tgt_EPSG, tgt_res, nodata_out, interp_method)

        # 5/ ELLIPSOID OR GEOID

        # raw data is ellipsoid WGS84
        if geoid:
            tmp_dem_geoid= os.path.join(tmp_dir, os.path.splitext(os.path.basename(tmp_dem))[0] + '_geoid.tif')
            ellipsoid_to_geoid(tmp_dem_proj,tmp_dem_geoid)
        else:
            tmp_dem_geoid=tmp_dem_proj

        list_tmp_dem_tomerge.append(tmp_dem_geoid)

    # 6/ MERGE ALL TILES

    tmp_dem_merged=os.path.join(tmp_dir,tile_name+'_merged.tif')
    if tag_merge:
        merge_rast_list(list_tmp_dem_tomerge,tmp_dem_merged)
    else:
        shutil.copytree(tmp_dir,os.path.join(os.path.dirname(dem_out),tile_name+'_subtiles'))

    # 7/ CLIP TO TILE EXTENT

    if not tag_clip:
        tmp_dem_clipped = os.path.join(tmp_dir,tile_name+'_clipped.tif')
        lat,lon= SRTMGL1_naming_to_latlon(tile_name)
        clip_rast_to_extent(tmp_dem_merged, tmp_dem_clipped, [lat, lon, lat + 1, lon + 1], 4326)
    else:
        tmp_dem_clipped = tmp_dem_merged

    shutil.move(tmp_dem_clipped,dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def REMA_mosaic_r1_1_tile(dir_REMA,tile_name,dem_out,filter_params=None,format_out='GTiff',tgt_EPSG=3031,tgt_res=None,nodata_out=-9999,interp_method=None,geoid=False,tag_lonlat_tile=False,path_tile_index=None,tag_merge=False,tag_clip=False):

    """
    :param dir_REMA: path to parent directory "8m" containing subdirectories of tar.gz archives (native FTP architecture)
    :param tile_name: either REMA tile name or 1x1° lat/lon tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param filter_params: filtering with REMA ERR file using rastlib.filter_nanarray function
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :param tag_lonlat_tile: True if tile_name follows SRTMGL1 tile naming, False if tile_name follows REMA tile naming
    :param path_tile_index: if tile_name is REMA format, specify path to ESRI REMA Tile Index
    :param tag_merge: if tile_name is REMA format, True to merge all ArcticDEM tiles to the 1x1° lat/lon extent
    :param tag_clip: if tile_name is REMA format, True to clip the 5x5° tile to the 1x1° lat/lon extent of tile_name
    :return:

    REMA release 1.1 product: ref:https://www.pgc.umn.edu/data/rema/

    Processing for 8m mosaic
    (100m, 500m and 1km versions are bundled in one .tif file)

    Tile name and processing is REMA tile naming convention by default
    Provide path to ESRI tile index file to use 1x1° lat/lon tiles and SRTMGL1 naming convention

    OPTIMAL DIRECTORY ARCHITECTURE: point to "8m" folder of similar architecture than: ftp://ftp.data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.0
    """

    # 1/ LOCATE TILE

    if not tag_lonlat_tile:
        subtile_dir=os.path.join(dir_REMA,tile_name)
        tile_tar_gz_list=[os.path.join(subtile_dir,tar_file) for tar_file in os.listdir(subtile_dir) if tar_file.endswith('.tar.gz')]
    else:
        lat_tile, lon_tile = SRTMGL1_naming_to_latlon(tile_name)
        extent = [lat_tile, lon_tile, lat_tile + 1, lon_tile + 1]
        # feature name in REMA_Tile_Index_Rel1.1
        feat_name = 'TILE'
        subtile_name_list=list_shp_field_inters_extent(path_tile_index, feat_name, extent,4326)
        subtile_dir_list = [os.path.join(dir_REMA,tile) for tile in subtile_name_list]
        tile_tar_gz_list=[]
        for i in range(len(subtile_dir_list)):
            tile_tar_gz_list=tile_tar_gz_list+[os.path.join(subtile_dir_list[i],tar_file) for tar_file in os.listdir(subtile_dir_list[i]) if tar_file.endswith('.tar.gz')]

    # 2/ EXTRACT TILE

    tmp_dir = create_tmp_dir_for_outfile(dem_out)

    list_tmp_dem = [os.path.join(tmp_dir, os.path.splitext(os.path.basename(tile_tar_gz))[0]+'_dem.tif') for tile_tar_gz in tile_tar_gz_list]
    for tile_tar_gz in tile_tar_gz_list:
        extract_file_from_tar_gz(tile_tar_gz,os.path.splitext(os.path.basename(tile_tar_gz))[0]+'_dem.tif',list_tmp_dem[tile_tar_gz_list.index(tile_tar_gz)])

    # list_tmp_err = [tmp_dir + os.path.splitext(os.path.basename(tile_tar_gz))[0]+'_err.tif' for tile_tar_gz in tile_tar_gz_list]
    for tile_tar_gz in tile_tar_gz_list:
        extract_file_from_tar_gz(tile_tar_gz,os.path.splitext(os.path.basename(tile_tar_gz))[0]+'_err.tif',list_tmp_dem[tile_tar_gz_list.index(tile_tar_gz)])

    list_tmp_dem_tomerge=[]
    for tmp_dem in list_tmp_dem:
        # 3/ FILTER TILE
        if filter_params is not None:

            tmp_err=tmp_dem[:-8]+'_err.tif'

            err = read_nanarray(tmp_err)
            _, filt = filter_nanarray(err, filter_params[0], filter_params[1], filter_params[2])

            dem = read_nanarray(tmp_dem)
            dem_filtered = np.array(dem)
            dem_filtered[filt] = np.NaN

            update_nanarray(tmp_dem, dem_filtered)

        # 4/ REPROJECT TILE

        # raw data is GeoTiff, 3031, 1 arc-sec and -9999 nodata_out
        if format_out == 'GTiff' and tgt_EPSG == 3031 and tgt_res is None and nodata_out is -9999:
            tmp_dem_proj = tmp_dem
        else:
            tmp_dem_proj = os.path.join(tmp_dir, os.path.splitext(os.path.basename(tmp_dem))[0] + '_proj.tif')
            warp_defaultUTM(tmp_dem, tmp_dem_proj, format_out, 3031, tgt_EPSG, tgt_res, nodata_out, interp_method)

        # 5/ ELLIPSOID OR GEOID

        # raw data is ellipsoid WGS84
        if geoid:
            tmp_dem_geoid= os.path.join(tmp_dir, os.path.splitext(os.path.basename(tmp_dem))[0] + '_geoid.tif')
            ellipsoid_to_geoid(tmp_dem_proj,tmp_dem_geoid)
        else:
            tmp_dem_geoid=tmp_dem_proj

        list_tmp_dem_tomerge.append(tmp_dem_geoid)

    # 6/ MERGE ALL TILES

    tmp_dem_merged=os.path.join(tmp_dir,tile_name+'_merged.tif')
    if tag_merge:
        merge_rast_list(list_tmp_dem_tomerge,tmp_dem_merged)
    else:
        shutil.copytree(tmp_dir,os.path.join(os.path.dirname(dem_out),tile_name+'_subtiles'))

    # 7/ CLIP TO TILE EXTENT

    if not tag_clip:
        tmp_dem_clipped = os.path.join(tmp_dir,tile_name+'_clipped.tif')
        lat,lon= SRTMGL1_naming_to_latlon(tile_name)
        clip_rast_to_extent(tmp_dem_merged, tmp_dem_clipped, [lat, lon, lat + 1, lon + 1], 4326)
    else:
        tmp_dem_clipped = tmp_dem_merged

    shutil.move(tmp_dem_clipped,dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def SRTMGL1_V003_tile(dir_SRTMGL1,tile_name,dem_out,filter_params=None,format_out='SRTMHGT',tgt_EPSG=4326,tgt_res=None,nodata_out=-32768,interp_method=None,geoid=True):

    """
    :param dir_SRTMGL1: path to directory directly containing zip of SRTMGL1 and SRTMGL1N
    :param tile_name: 1x1° tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param filter_params: filtering with SRTMGL1 NUM file using rastlib.filter_nanarray function
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :return:

    SRTMGL1: Shuttle Radar Topography Mission Global 1-arcsec V003 from NASA at: https://lpdaac.usgs.gov/sites/default/files/public/measures/docs/NASA_SRTM_V3.pdf

    NUM file description
    1 = Water-masked SRTM void *
    2 = Water-masked SRTM non-void *
    5 = GDEM elevation = 0 in SRTM void (helped correct ocean masking)
    11 = NGA-interpolated SRTM (were very small voids in SRTMv1) 21 = GMTED2010 oversampled from 7.5 arc-second postings
    25 = SRTM within GDEM **
    31 = NGA fill of SRTM via GDEM***
    51 = USGS NED
    52 = USGS NED via GDEM
    53 = Alaska USGS NED via GDEM
    72 = Canadian Digital Elevation Data (CDED) via GDEM 101-200 = ASTER scene count (count limited to 100)
    201-224 = SRTM swath count (non-voided swaths) Actual maximum = 24

    OPTIMAL DIRECTORY ARCHITECTURE: SRTMGL1 and SRTMGL1N zip archives in the same directory
    """

    #1/ LOCATE 1x1° tile

    #nomenclature of SRTMGL1 and SRTMGL1N zip files
    hgt_zip=os.path.join(dir_SRTMGL1,tile_name+'.SRTMGL1.hgt.zip')
    num_zip=os.path.join(dir_SRTMGL1,tile_name+'.SRTMGL1N.num.zip') #NOTE: nomenclature of V002 was SRTMGL1.num.zip for the num file

    #2/ EXTRACT TILE

    tmp_dir=create_tmp_dir_for_outfile(dem_out)

    zip_ref = zipfile.ZipFile(hgt_zip, 'r')
    zip_ref.extractall(tmp_dir)
    zip_ref.close()

    zip_ref = zipfile.ZipFile(num_zip, 'r')
    zip_ref.extractall(tmp_dir)
    zip_ref.close()

    tmp_hgt=os.path.join(tmp_dir,tile_name+'.hgt')
    tmp_num=os.path.join(tmp_dir,tile_name+'.num')

    #3/ FILTER TILE
    if filter_params is not None:

        num = read_nanarray(tmp_num)
        _, filt = filter_nanarray(num,filter_params[0],filter_params[1],filter_params[2])

        hgt = read_nanarray(tmp_hgt)
        hgt_filtered = np.array(hgt)
        hgt_filtered[filt]=np.NaN

        update_nanarray(tmp_hgt,hgt_filtered)

    #4/ REPROJECT TILE

    # raw data is SRTMHGT, 4326, 1 arc-sec and -32768 nodata_out
    if format_out == 'SRTMHGT' and tgt_EPSG == 4326 and tgt_res is None and nodata_out is -32768:
        tmp_hgt_proj = tmp_hgt
    else:
        tmp_hgt_proj = os.path.join(tmp_dir, tile_name + '_proj.tif')
        warp_defaultUTM(tmp_hgt, tmp_hgt_proj, format_out, 4326, tgt_EPSG, tgt_res, nodata_out, interp_method)

    #5/ ELLIPSOID OR GEOID

    # raw data is geoid EGM96
    if not geoid:
        geoid_to_ellipsoid(tmp_hgt_proj, dem_out)
    else:
        shutil.move(tmp_hgt_proj,dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def ASTGTM2_V002_tile(dir_ASTGTM2,tile_name,dem_out,filter_params=None,format_out='GTiff',tgt_EPSG=4326,tgt_res=None,nodata_out=-9999,interp_method=None,geoid=True):


    """
    :param dir_ASTGTM2: path to directory directly containing zip of ASTGTM2
    :param tile_name: 1x1° tile name (SRTMGL1/classic naming convention)
    :param dem_out: path to DEM out file
    :param filter_params: filtering with ASTGTM2 NUM file using rastlib.filter_nanarray function
    :param format_out: output format, GDAL naming (e.g.: 'GTiff','HDF4', ...) ; see: https://www.gdal.org/formats_list.html
    :param tgt_EPSG: EPSG of output projection
    :param tgt_res: output resolution, GDAL naming [xres, yres]
    :param nodata_out: output no-data value
    :param interp_method: resampling method, GDAL naming 'bilinear', 'neir', 'cubic', etc..
    :param geoid: True, converts to geoid if is ellipsoid; False converts to ellipsoid if is geoid
    :return:

    ASTGTM2: ASTER Global Digital Elevation Model V002 from NASA at: https://lpdaac.usgs.gov/node/1079

    NUM file description
    Fill value: -1 to -11
    -1	SRTM3 V003
    -2	SRTM3 V002
    -5	NED
    -6	CDED
    -11	Alaska DEM
    Valid range (ASTER count):
    0 to 200

    OPTIMAL DIRECTORY ARCHITECTURE: ASTGTM2 zip archives in the same directory
    """

    #1/ LOCATE 1x1° tile

    #nomenclature of ASTGTM2 zip files
    tile_zip=os.path.join(dir_ASTGTM2,'ASTGTM2_'+tile_name+'.zip')

    dem_file = 'ASTGTM2_'+tile_name+'_dem.tif'
    num_file = 'ASTGTM2_'+tile_name+'_num.tif'

    #2/ EXTRACT TILE

    tmp_dir=create_tmp_dir_for_outfile(dem_out)

    tmp_dem=os.path.join(tmp_dir,dem_file)
    tmp_num=os.path.join(tmp_dir,num_file)

    extract_file_from_zip(tile_zip,dem_file,tmp_dem)
    extract_file_from_zip(tile_zip,num_file,tmp_num)

    #3/ FILTER TILE
    if filter_params is not None:

        num = read_nanarray(tmp_num)
        _, filt = filter_nanarray(num,filter_params[0],filter_params[1],filter_params[2])

        dem = read_nanarray(tmp_dem)
        dem_filtered = np.array(dem)
        dem_filtered[filt]=np.NaN

        update_nanarray(tmp_dem,dem_filtered)

    #4/ REPROJECT TILE

    # raw data is GTiff, 4326, 1 arc-sec and -32768 nodata_out
    if format_out == 'GTiff' and tgt_EPSG == 4326 and tgt_res is None and nodata_out is -9999:
        tmp_dem_proj = tmp_dem
    else:
        tmp_dem_proj = os.path.join(tmp_dir, tile_name + '_proj.tif')
        warp_defaultUTM(tmp_dem, tmp_dem_proj, format_out, 4326, tgt_EPSG, tgt_res, nodata_out, interp_method)

    #5/ ELLIPSOID OR GEOID

    # raw data is geoid EGM96
    if not geoid:
        geoid_to_ellipsoid(tmp_dem_proj,dem_out)
    else:
        shutil.move(tmp_dem_proj,dem_out)

    remove_tmp_dir_for_outfile(dem_out)


def NASADEM_tile():

    """
    Waiting for finalized product with header...
    """