# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

TILE DIVISION LIBRARY:
Library of Python functions for manipulating vectors/rasters/files for 1x1° tile based calculations/file manipulation

LIST:

- copy/rename/move/delete subfolders

- copy/rename/move/delete files within subfolders

- update files within subfolders

- count files within subfolders
"""
from __future__ import print_function
from osgeo import gdal, ogr, osr
# from tqdm import tqdm
import numpy as np
from shlib import create_tmp_dir_for_outfile, remove_tmp_dir_for_outfile
from misclib import latlon_to_UTM, SRTMGL1_naming_to_latlon, extended_factor_latitude_L1A, along_track_Terra
from rastlib import extent_rast, poly_from_extent, merge_rast_list, read_proj_rast
from vectlib import proj_WKT_to_EPSG, inters_list_poly_with_poly, coord_trans_wkt_or_EPSG, get_poly_centroid, extent_from_poly

def extend_tile_raster_list(list_raster_in,extend_factor,extend_epsg,list_raster_out,nodata_val):

    list_poly_extended=[]
    list_poly=[]
    list_epsg=[]

    for raster_in in list_raster_in:
        tmp_extent,tmp_proj_wkt=extent_rast(raster_in)
        tmp_poly = poly_from_extent(tmp_extent)

        #detect UTM EPSG of central raster
        transform = coord_trans_wkt_or_EPSG(True, tmp_proj_wkt, False, extend_epsg)
        tmp_poly.Transform(transform)
        center_lon, center_lat = get_poly_centroid(tmp_poly)
        tmp_epsg, _ = latlon_to_UTM(center_lat, center_lon)
        tmp_extent_proj = extent_from_poly(tmp_poly)

        list_poly.append(tmp_poly)
        list_epsg.append(int(tmp_epsg))

        tmp_extended_extent=[tmp_extent_proj[0]-extend_factor[0],tmp_extent_proj[1]-extend_factor[1],tmp_extent_proj[2]+extend_factor[0],tmp_extent_proj[3]+extend_factor[1]]
        list_poly_extended.append(poly_from_extent(tmp_extended_extent))

    list_list_inters_extended=[]
    for poly in list_poly_extended:
        list_inters=inters_list_poly_with_poly(list_poly,poly)
        list_list_inters_extended.append(list_inters)

    for raster_in in list_raster_in:

        ind=list_raster_in.index(raster_in)
        list_inters_extended = list_list_inters_extended[ind]

        print('FILE NUMBER ' + str(ind) + ': '+list_raster_out[ind])
        list_inters_raster = [list_raster_in[list_poly.index(inters_extended)] for inters_extended in list_inters_extended]
        print(list_inters_raster)
        merge_rast_list(list_inters_raster,list_raster_out[ind],tgt_EPSG=list_epsg[ind],nodata_in=nodata_val,nodata_out=nodata_val)


def area_intersect_geom_listpoly(list_tiles,list_poly,geom,type_geom):

    #WARNING: reprojected polygons are updated in the original list

    list_tiles_inters=[]
    list_area_inters=[]
    list_tot_inters=[]

    for i in range(len(list_tiles)):

        # progress_bar(float(i+1),float(len(list_tiles)),20)

        if type_geom=='Polygon':

            tmpintersect = geom.Intersection(list_poly[i])  #polygon intersection with extent

            if not tmpintersect.IsEmpty():

                list_tiles_inters.append(list_tiles[i]) #add tile to list if intersection is not empty

        elif type_geom=='GeometryCollection':

            # reproj multi-polygon
            lat, lon = SRTMGL1_naming_to_latlon(list_tiles[i])
            epsg, utm_zone = latlon_to_UTM(lat + 0.5, lon + 0.5)

            source = osr.SpatialReference()
            source.ImportFromEPSG(4326)

            target = osr.SpatialReference()
            target.ImportFromEPSG(int(epsg))

            transform = osr.CoordinateTransformation(source, target)
            tile_area = 0

            # tmpclip = extent_geom.Clip(list_poly)

            for g in geom:

                tmpintersect = g.Intersection(list_poly[i])  #polygon intersection with extent

                if not tmpintersect.IsEmpty():
                    tmpintersect.Transform(transform)
                    tile_area = tile_area + tmpintersect.GetArea()

            if tile_area > 0:
                tmppoly = list_poly[i]
                tmppoly.Transform(transform)
                tot_area = tmppoly.GetArea()
                list_tiles_inters.append(list_tiles[i])  #add tile to list if intersection is not empty
                list_area_inters.append(tile_area)
                list_tot_inters.append(tot_area)

    return list_tiles_inters,list_area_inters,list_tot_inters


def stack_tile_polygon(list_tiles,flag_extended):

    if list_tiles is None:

        print('Creating all 1x1 degree polygons globally...')

        list_tiles=[]

        #naming convention of SRTM: longitude=W180 to W1 then E000 to E179 ; latitude=S1 to S90 then N00 to N89
        for ymin in np.arange(-90,90): #iterating from -90 to 89
            for xmin in np.arange(-180,180): #iterating from -180 to 179
                if xmin<0:
                    if ymin<0:
                        list_tiles.append('S'+str(int(np.abs(np.round(ymin)))).zfill(2)+"W"+str(int(np.abs(np.round(xmin)))).zfill(3))
                    else:
                        list_tiles.append('N'+str(int(np.abs(np.round(ymin)))).zfill(2)+"W"+str(int(np.abs(np.round(xmin)))).zfill(3))
                else:
                    if ymin<0:
                        list_tiles.append('S'+str(int(np.abs(np.round(ymin)))).zfill(2)+"E"+str(int(np.abs(np.round(xmin)))).zfill(3))
                    else:
                        list_tiles.append('N'+str(int(np.abs(np.round(ymin)))).zfill(2)+"E"+str(int(np.abs(np.round(xmin)))).zfill(3))

    #stacks coordinates of corresponding tile polygon in 3dim array

    nb_tiles = len(list_tiles)

    coords = np.zeros([nb_tiles, 5, 2])

    for i in range(nb_tiles):

        lat, lon = SRTMGL1_naming_to_latlon(list_tiles[i])
        lat=lat+0.5
        lon=lon+0.5

        if flag_extended:
            ext_lon, ext_lat = extended_factor_latitude_L1A(lat,abs(along_track_Terra(lat)))
        else:
            # extension from center of tile
            ext_lon = 0.5
            ext_lat = 0.5

        # SW corner
        coords[i, 0, 0] = lon + ext_lon
        coords[i, 0, 1] = lat + ext_lat
        # NW corner
        coords[i, 1, 0] = lon - ext_lon
        coords[i, 1, 1] = lat + ext_lat
        # NE corner
        coords[i, 2, 0] = lon - ext_lon
        coords[i, 2, 1] = lat - ext_lat
        # SE corner
        coords[i, 3, 0] = lon + ext_lon
        coords[i, 3, 1] = lat - ext_lat
        # closing polygon: SW corner
        coords[i, 4, 0] = lon + ext_lon
        coords[i, 4, 1] = lat + ext_lat

    print('Creating and stacking 1x1 degree polygons...')
    #reading coord array shape
    [a, _, __] = np.shape(coords)

    list_poly=[]

    for i in range(a):
        listcoord = zip(coords[i, :, 0].tolist(), coords[i, :, 1].tolist()) #zip coordinate array in list of tuples

        ring = ogr.Geometry(ogr.wkbLinearRing) #creating polygon ring
        for coord in listcoord:
            ring.AddPoint(coord[0], coord[1])

        tmppoly = ogr.Geometry(ogr.wkbPolygon)
        tmppoly.AddGeometry(ring) #creating polygon

        list_poly.append(tmppoly) #stacking polygon as a list

    return list_tiles, list_poly