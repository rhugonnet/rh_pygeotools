# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

Python functions for vector and geometry manipulation with OGR/OSR Python bindings

LIST OF FUNCTIONS:

- read_geom_from_shp: read geometry from ESRI shapefile

- write_geom_to_shp: write (polygon or multipolygon) geometry to ESRI shapefile

- create_valid_multipoly_from_shp: create multipolygon and correct the geometries with a "0" buffer from an ESRI shapefile

- extent_shp: extract extent from an ESRI shapefile

- area_intersect_geom_listpoly: calculate intersecting area of a geometry (polygon or geomcollection) and a polygon stack
# !warning: input polygon list is reprojected during the calculation

- polygon_list_to_multipoly: create multipolygon from polygon list

- union_cascaded_multipoly: derive cascaded union of multipolygon

- intersect_poly: derive intersection of polygon geometries

- coord_trans_wkt_or_EPSG: create a transform object to reproject geometries based on source/target projections, imported from either EPSG or WKT (format read by gdal/ogr using GetProjection())

- clip_shp_to_extent:

- merge_shp_list: merge all ESRI shapefiles (.shp) in a list

"""
from __future__ import print_function
import os, shutil, sys
from subprocess import Popen
from osgeo import ogr, osr, gdal
from misclib import latlon_to_UTM
from shapely import geometry
from shapely.algorithms import polylabel
import json

def extent_shp_ref(fn_shp_in,fn_ref):

    ext, proj = extent_shp(fn_shp_in)
    poly = poly_from_extent(ext)

    ds = gdal.Open(fn_ref,gdal.GA_ReadOnly)
    proj_ref = ds.GetProjection()

    trans =coord_trans_wkt_or_EPSG(True,proj.ExportToWkt(),True,proj_ref)

    poly.Transform(trans)

    ext= extent_from_poly(poly)

    return ext

def copy_shp_fn(fn_shp_in,fn_shp_out):

    os.system('ogr2ogr '+fn_shp_out+' '+fn_shp_in)

def isempty_firstfeat(fn_shp):

    ds_shp = ogr.Open(fn_shp)
    layer = ds_shp.GetLayer(0)
    feature = layer.GetFeature(0)
    if feature is None:
        tag = True
    else:
        poly = feature.GetGeometryRef()
        if poly is None:
            tag = True
        else:
            tag = False

    return tag

def poi_polygon(fn_shp,fn_poi_shp,tolerance):

    ds_shp = ogr.Open(fn_shp)
    layer = ds_shp.GetLayer(0)
    proj_wkt=layer.GetSpatialRef().ExportToWkt()
    feature = layer.GetFeature(0)
    export = feature.GetGeometryRef().ExportToJson()

    export_valid = json.loads(export)
    geom = geometry.shape(export_valid)

    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)

    poi = polylabel.polylabel(geom,tolerance)

    write_poly_to_shp(poi.wkt,fn_poi_shp,srs=srs,layer_type=ogr.wkbPoint,export_wkt=False)

    ds_shp = None


def clip_shp_to_shp(fn_shp_toclip,fn_shp_clipping,fn_shp_out):

    os.system('ogr2ogr -clipsrc ' + fn_shp_clipping + ' -f "ESRI Shapefile" ' + fn_shp_out + ' ' + fn_shp_toclip)

def extract_list_feat_ogr(fn_shp_in,fn_shp_out,feat_name,list_feat_id):

    layer_name=os.path.splitext(os.path.basename(fn_shp_in))[0]
    os.system('ogr2ogr -sql "'+"SELECT * FROM "+layer_name+" WHERE "+feat_name+" IN ('"+"','".join(list_feat_id)+"'"+')" '+ fn_shp_out+' '+fn_shp_in)

def extract_feat_as_shp(fn_shp_in,fn_shp_out,feat_name,feat_id):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds_shp = driver.Open(fn_shp_in, 0)
    layer = ds_shp.GetLayer()
    proj_shp = layer.GetSpatialRef().ExportToWkt()
    geom_out=None
    for feature in layer:
        geom = feature.GetGeometryRef()
        feat_val = feature.GetField(feat_name)
        if feat_val==feat_id:
            geom_out = geom
            break

    if geom_out is None:
        print('Could not find feature id.')
        sys.exit()

    srs=osr.SpatialReference()
    srs.ImportFromWkt(proj_shp)
    write_poly_to_shp(geom_out,fn_shp_out,srs=srs)

def union_shp_fn(fn_shp_1,fn_shp_2,fn_shp_out):

    geom1, proj_wkt = read_geom_from_shp(fn_shp_1)
    geom2, proj_wkt2 = read_geom_from_shp(fn_shp_2)

    # poly1 = geomcol_to_valid_multipoly(geom1)
    # poly2 = geomcol_to_valid_multipoly(geom2)

    transform = coord_trans_wkt_or_EPSG(True, proj_wkt2, True, proj_wkt)
    geom2.Transform(transform)

    inters = geom1.Union(geom2)

    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)
    write_poly_to_shp(inters,fn_shp_out,srs=srs)

def inters_shp_fn(fn_shp_1,fn_shp_2,fn_shp_out):

    geom1, proj_wkt = read_geom_from_shp(fn_shp_1)
    geom2, proj_wkt2 = read_geom_from_shp(fn_shp_2)

    transform = coord_trans_wkt_or_EPSG(True, proj_wkt2, True, proj_wkt)

    geom2.Transform(transform)
    geom1 = ogr.ForceToMultiPolygon(geom1)
    geom2 = ogr.ForceToMultiPolygon(geom2)
    inters = geom1.Intersection(geom2)

    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)
    write_poly_to_shp(inters,fn_shp_out,srs=srs)

def buffer_shp_fn(fn_shp_in,fn_shp_out,buffer,meters=True):

    geom, proj_wkt = read_geom_from_shp(fn_shp_in)
    if meters:
        #project to UTM of centroid by default
        transform = coord_trans_wkt_or_EPSG(True, proj_wkt, False, 4326)
        geom.Transform(transform)
        center_lon, center_lat = get_poly_centroid(geom)
        epsg, utm_zone = latlon_to_UTM(center_lat, center_lon)
        transform2 = coord_trans_wkt_or_EPSG(False, 4326, False, int(epsg))
        geom.Transform(transform2)
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(int(epsg))
    else:
        srs = osr.SpatialReference()
        srs.ImportFromWkt(proj_wkt)

    buff = geom.Buffer(buffer)

    write_poly_to_shp(buff,fn_shp_out,srs=srs)

def simplify_shp_fn(fn_shp_in,fn_shp_out,tolerance):

    os.system('ogr2ogr -simplify '+str(int(tolerance))+' '+fn_shp_out+' '+fn_shp_in)

    #or:

    # geomcol,proj_wkt = read_geom_from_shp(fn_shp_in)
    # geomcol = ogr.ForceToMultiPolygon(geomcol)
    # simplified = geomcol.SimplifyPreserveTopology(tolerance)
    #
    # srs = osr.SpatialReference()
    # srs.ImportFromWkt(proj_wkt)
    # write_poly_to_shp(simplified,fn_shp_out,srs=srs,layer_type=ogr.wkbMultiPolygon)

def inters_list_poly_with_poly(list_poly,poly):

    list_inters=[]
    for poly_2 in list_poly:
        inters = poly_2.Intersection(poly)

        if not inters.IsEmpty():
            list_inters.append(poly_2)

    return list_inters

def proj_WKT_to_EPSG(proj_wkt):

    proj=osr.SpatialReference()
    proj.ImportFromWkt(proj_wkt)
    epsg = proj.AutoIdentifyEPSG()

    return epsg

def list_shp_field_inters_extent(in_shp,field_name,extent,extent_EPSG):

    poly = poly_from_extent(extent)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(in_shp, 0)
    layer = ds.GetLayer()

    proj_shp = layer.GetSpatialRef().ExportToWkt()

    trans = coord_trans_wkt_or_EPSG(False,extent_EPSG,True,proj_shp)

    poly.Transform(trans)

    list_field_inters=[]
    for feat in layer:
        feat_geom = feat.GetGeometryRef()
        inters = feat_geom.Intersection(poly)

        if not inters.IsEmpty():
            list_field_inters.append(feat.GetField(field_name))

    return list_field_inters

def get_poly_centroid(poly):

    centroid = poly.Centroid()

    center_lon, center_lat, _ = centroid.GetPoint()

    return center_lon, center_lat

def read_proj_shp(in_shp):
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataset = driver.Open(in_shp)

    #from layer
    layer = dataset.GetLayer()
    proj_wkt = layer.GetSpatialRef()

    return proj_wkt.ExportToWkt()

def merge_shp_list(list_shp,out_merged):

    tmp_dir=os.path.join(os.path.dirname(out_merged),'tmp_dir')+os.path.sep

    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    for i in range(len(list_shp)-1):
        if i==0:
            in_1=list_shp[i]
        else:
            in_1=out_shp

        in_2=list_shp[i+1]

        out_shp=tmp_dir+'tmp_'+str(i)+'.shp'

        os.system('ogr2ogr -f "ESRI Shapefile" '+out_shp+' '+in_1)
        os.system('ogr2ogr -f "ESRI Shapefile" -append -update '+out_shp+' '+in_2)

    os.system('ogr2ogr '+out_merged+' '+out_shp)

    shutil.rmtree(tmp_dir)


def clip_shp_to_extent(in_shp,ext_shp,out_shp):

    source_ds = ogr.Open(ext_shp)
    source_layer = source_ds.GetLayer()
    xmin, xmax, ymin, ymax = source_layer.GetExtent()

    os.system('ogr2ogr -clipsrc '+str(xmin)+' '+str(ymin)+' '+str(xmax)+' '+str(ymax)+' '+out_shp+' '+in_shp)


def compare_proj_wkt_or_EPSG(tag_src_wkt,proj_src,tag_tgt_wkt,proj_tgt):

    # choice between most used projection formats: Wkt or EPSG
    source_proj = osr.SpatialReference()
    if tag_src_wkt:
        source_proj.ImportFromWkt(proj_src)
        source_epsg=source_proj.AutoIdentifyEPSG()
    else:
        source_epsg=proj_src

    target_proj = osr.SpatialReference()
    if tag_tgt_wkt:
        target_proj.ImportFromWkt(proj_tgt)
        target_epsg=target_proj.AutoIdentifyEPSG()
    else:
        target_epsg=proj_tgt

    return source_epsg == target_epsg, source_epsg, target_epsg


def coord_trans_wkt_or_EPSG(tag_src_wkt,proj_src,tag_tgt_wkt,proj_tgt):

    #choice between most used projection formats: Wkt or EPSG
    source_proj = osr.SpatialReference()
    if tag_src_wkt:
        source_proj.ImportFromWkt(proj_src)
    else:
        source_proj.ImportFromEPSG(proj_src)

    target_proj = osr.SpatialReference()
    if tag_tgt_wkt:
        target_proj.ImportFromWkt(proj_tgt)
    else:
        target_proj.ImportFromEPSG(proj_tgt)

    transform = osr.CoordinateTransformation(source_proj, target_proj)

    return transform


def intersect_poly(poly1,poly2):

    inters=poly1.Intersection(poly2)

    return inters


def union_cascaded_multipoly(multipoly):

    print('Calculating cascaded union of multipolygon...')
    # calculating cascaded union of multipolygon
    cascadedpoly = multipoly.UnionCascaded()

    return cascadedpoly


def polygon_list_to_multipoly(list_poly):

    print('Creating multipolygon from polygon list...')

    #create empty multipolygon
    multipoly = ogr.Geometry(ogr.wkbMultiPolygon)

    for i in range(len(list_poly)):
        #stacking polygons in multipolygon
        multipoly.AddGeometry(list_poly[i])

    return multipoly

def geom_extent(geom):

    env = geom.GetEnvelope()
    extent = [env[0], env[2], env[1], env[3]]

    return extent

def poly_from_extent(extent):

    xmin, ymin, xmax, ymax = extent

    ring = ogr.Geometry(ogr.wkbLinearRing)  # creating polygon ring
    ring.AddPoint(xmin, ymin)
    ring.AddPoint(xmax, ymin)
    ring.AddPoint(xmax, ymax)
    ring.AddPoint(xmin, ymax)
    ring.AddPoint(xmin, ymin)

    poly = ogr.Geometry(ogr.wkbPolygon)
    poly.AddGeometry(ring)  # creating polygon

    return poly

def extent_from_poly(poly):
    env = poly.GetEnvelope()

    extent = env[0], env[2], env[1], env[3]

    return extent

def extent_shp(in_shp):

    #get layer from ESRI shapefile
    inShapefile = in_shp
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(inShapefile, 0)
    inLayer = inDataSource.GetLayer()
    ext = inLayer.GetExtent()
    proj_wkt = inLayer.GetSpatialRef()

    xmin=min(ext[0],ext[1])
    ymin=min(ext[2],ext[3])
    xmax=max(ext[0], ext[1])
    ymax=max(ext[2], ext[3])

    extent = [xmin, ymin, xmax, ymax]

    return extent, proj_wkt

def convhull_of_geomcol(geomcol):

    #calculate convex hull of collection
    convhull = geomcol.ConvexHull()

    return convhull


def geomcol_to_valid_multipoly(geomcol,flag_buffer=False):

    multi = ogr.Geometry(ogr.wkbMultiPolygon)
    for g in geomcol:
        multi.AddGeometry(g)

    if flag_buffer:
        print('Correcting shapefile to consider only valid geometries...')
        #removing invalid geometries of multi-polygons (intersecting ring points) by using a '0' buffer
        multi_valid=multi.Buffer(0)

        return multi_valid
    else:
        return multi


def read_geom_from_shp(in_shp):

    print('Reading shapefile: ' + in_shp + '...')

    #get layer from ESRI shapefile
    inShapefile = in_shp
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(inShapefile, 0)
    inLayer = inDataSource.GetLayer()
    proj = inLayer.GetSpatialRef()


    #collect all geometry in a geometry collection
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in inLayer:
        geomcol.AddGeometry(feature.GetGeometryRef())

    return geomcol, proj.ExportToWkt()


def write_poly_to_shp(geom,out_shp,srs=None,layer_name='NA',field_id='ID',field_val='123',layer_type=ogr.wkbPolygon,export_wkt=True,flag_zip=False):

    if os.path.exists(out_shp):
        shutil.rmtree(out_shp,ignore_errors=True)

    print('Writing shapefile to disk: ' + out_shp +'...')

    """
    https://gis.stackexchange.com/a/52708/8104
    """
    if srs is None:
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

    #convert to a shapefile with OGR
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(out_shp)
    layer = ds.CreateLayer(layer_name, srs, layer_type)
    # Add one attribute
    layer.CreateField(ogr.FieldDefn(field_id ,ogr.OFTInteger))
    defn = layer.GetLayerDefn()

    #if there are multiple geometries, put the "for" loop here
    #create a new feature (attribute and geometry)
    feat = ogr.Feature(defn)
    feat.SetField(field_id, field_val)

    #make geometry
    if export_wkt:
        geom = ogr.CreateGeometryFromWkt(geom.ExportToWkt())
    else:
        geom = ogr.CreateGeometryFromWkt(geom)
    feat.SetGeometry(geom)

    layer.CreateFeature(feat)
    feat = geom = None

    ds = layer = feat = geom = None

    if flag_zip:
        if out_shp[-1:]=='/':
            out_shp=out_shp[:-1]

        shutil.make_archive(os.path.abspath(os.path.join(out_shp,os.pardir))+'/'+os.path.basename(out_shp),'zip',out_shp)
