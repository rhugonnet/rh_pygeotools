# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

RASTER LIBRARY:
Library of Python functions for raster manipulation

TODO: gdal.Rasterize Python bindings have been corrected, should now be able to specify the type properly: need to change everything to boolean when masking

-
"""
from __future__ import print_function
import os, sys
import numpy as np
from osgeo import gdal, gdalconst, ogr, osr
from misclib import latlon_to_UTM
from vectlib import coord_trans_wkt_or_EPSG, compare_proj_wkt_or_EPSG, read_proj_shp, poly_from_extent, get_poly_centroid, extent_from_poly
from shlib import create_tmp_dir_for_outfile, remove_tmp_dir_for_outfile
import shutil
import operator
gdal.UseExceptions()


def list_valid_feat_intersect(fn_raster, fn_shp, feat_name, perc_min_pixel_valid=80.):

    # first, get intersection list of features with raster extent
    extent,proj_wkt=extent_rast(fn_raster)
    poly=poly_from_extent(extent)

    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds_shp = driver.Open(fn_shp, 0)
    layer = ds_shp.GetLayer()
    proj_shp = layer.GetSpatialRef().ExportToWkt()

    trans = coord_trans_wkt_or_EPSG(True,proj_wkt,True,proj_shp)
    poly.Transform(trans)

    list_feat_inters=[]
    for feature in layer:
        geom = feature.GetGeometryRef()
        inters= geom.Intersection(poly)
        if not inters.IsEmpty():
            feat_val = feature.GetField(feat_name)
            list_feat_inters.append(feat_val)
    ds_shp.Destroy()

    #now, rasterize those features and find out valid pixel intersecting
    #open shp with GDAL this time
    ds_shp = gdal.OpenEx(fn_shp, gdal.OF_VECTOR)
    rast = read_nanarray(fn_raster)
    layer_name = os.path.splitext(os.path.basename(fn_shp))[0]
    list_feat_valid=[]
    # list_perc_inters=[]
    for feat_inters in list_feat_inters:
        ds_raster = create_mem_raster_on_ref(fn_raster)
        rasterize_feat_shp_ds(ds_shp,layer_name,feat_name,feat_inters,ds_raster)
        mask = (ds_raster.GetRasterBand(1).ReadAsArray() ==1)

        feat_rast = rast[mask]
        nb_pixel_tot=len(feat_rast)
        nb_pixel_valid=len(feat_rast[np.isfinite(feat_rast)])

        if nb_pixel_valid>perc_min_pixel_valid/100.*nb_pixel_tot:
            list_feat_valid.append(feat_inters)

        # list_feat_valid.append(feat_inters)
        # list_perc_inters.append(float(nb_pixel_valid/nb_pixel_tot))
        ds_raster=None

    ds_shp=None

    return list_feat_valid

def polygonize_fn(fn_raster_in,fn_shp_out):

    src_ds = gdal.Open(fn_raster_in,gdalconst.GA_ReadOnly)
    src_band = src_ds.GetRasterBand(1)

    mask_band = src_band.GetMaskBand()

    srs=osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjectionRef())

    #remove the .shp if specified by user
    # if fn_shp_out[-4:] == '.shp':
    #     fn_shp_out = fn_shp_out[:-4]

    drv = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(fn_shp_out):
        drv.DeleteDataSource(fn_shp_out)
    dst_ds=drv.CreateDataSource(fn_shp_out)
    dst_layer= dst_ds.CreateLayer('polygonized',srs=srs)
    # newField = ogr.FieldDefn('GLWD_all', ogr.OFTInteger)
    # dst_layer.CreateField(newField)

    gdal.Polygonize(src_band,mask_band,dst_layer,-1,['8CONNECTED=8'],callback=None)
    dst_ds.Destroy()
    src_ds = None

def proximity_rast_fn(fn_raster_in,fn_raster_out,val=1):

    ds_raster_in=gdal.Open(fn_raster_in,gdalconst.GA_ReadOnly)

    ds_proxy = proximity_rast_ds(ds_raster_in,val=val)

    write_nanarray(fn_raster_out,fn_raster_in,ds_proxy.GetRasterBand(1).ReadAsArray())

def proximity_rast_ds(ds_raster_in,val=1):

    drv = gdal.GetDriverByName('MEM')
    proxy_ds = drv.Create('', ds_raster_in.RasterXSize, ds_raster_in.RasterYSize, 1, gdal.GetDataTypeByName('Float32'))

    gdal.ComputeProximity(ds_raster_in.GetRasterBand(1), proxy_ds.GetRasterBand(1), ["VALUES="+str(val), "DISTUNITS=GEO"])

    return proxy_ds

def proximity_shp(fn_shp,raster_ref,type_prox='both'):

    # create memory raster based on reference raster
    raster_ds = create_mem_raster_on_ref(raster_ref)
    # open .shp with GDAL
    shp_ds = gdal.OpenEx(fn_shp, gdal.OF_VECTOR)
    # rasterize
    opts = gdal.RasterizeOptions(burnValues=[1], bands=[1])
    gdal.Rasterize(raster_ds, shp_ds, options=opts)
    # close shp
    shp_ds = None

    drv = gdal.GetDriverByName('MEM')
    proxy_ds = drv.Create('', raster_ds.RasterXSize, raster_ds.RasterYSize, 1, gdal.GetDataTypeByName('Float32'))

    if type_prox=='exterior':

        gdal.ComputeProximity(raster_ds.GetRasterBand(1),proxy_ds.GetRasterBand(1),["VALUES=1","DISTUNITS=GEO"])

        proxy_array= proxy_ds.GetRasterBand(1).ReadAsArray()

    elif type_prox=='interior':

        raster_arr = raster_ds.GetRasterBand(1).ReadAsArray()
        mask = (raster_arr == 1)
        raster_arr[mask] = raster_ds.GetRasterBand(1).GetNoDataValue()
        raster_arr[~mask] = 1
        raster_ds.GetRasterBand(1).WriteArray(raster_arr)

        gdal.ComputeProximity(raster_ds.GetRasterBand(1),proxy_ds.GetRasterBand(1),["VALUES=1","DISTUNITS=GEO"])

        proxy_array= proxy_ds.GetRasterBand(1).ReadAsArray()

    elif type_prox=='both':

        gdal.ComputeProximity(raster_ds.GetRasterBand(1),proxy_ds.GetRasterBand(1),["VALUES=1","DISTUNITS=GEO"])

        proxy_ext = proxy_ds.GetRasterBand(1).ReadAsArray()

        raster_arr = raster_ds.GetRasterBand(1).ReadAsArray()
        mask = (raster_arr == 1)
        raster_arr[mask] = raster_ds.GetRasterBand(1).GetNoDataValue()
        raster_arr[~mask] = 1
        raster_ds.GetRasterBand(1).WriteArray(raster_arr)

        gdal.ComputeProximity(raster_ds.GetRasterBand(1), proxy_ds.GetRasterBand(1), ["VALUES=1", "DISTUNITS=GEO"])

        proxy_int = proxy_ds.GetRasterBand(1).ReadAsArray()

        proxy_array = proxy_ext + proxy_int

    else:
        sys.exit('Type of proximity:'+ type_prox+ ' not recognized. Must be "interior", "exterior" or "both".')

    return proxy_array

def rasterize_feat_shp_fn(fn_shp,fn_raster_ref,fn_mask_out,feat_name,feat_id,exclude=False,all_touched=False):

    raster_ds = create_mem_raster_on_ref(fn_raster_ref)
    shp_ds = gdal.OpenEx(fn_shp, gdal.OF_VECTOR)

    rasterize_feat_shp_ds(shp_ds, os.path.splitext(os.path.basename(fn_shp))[0], feat_name, feat_id, raster_ds,exclude=exclude,all_touched=all_touched)

    mask_out= (raster_ds.GetRasterBand(1).ReadAsArray() == 1)
    write_nanarray(fn_mask_out,fn_raster_ref,mask_out)
    raster_ds = None
    shp_ds=None


def rasterize_shp(fn_shp,fn_raster_ref,all_touched=False):

    #create memory raster based on reference raster
    raster_ds = create_mem_raster_on_ref(fn_raster_ref)
    #open .shp with GDAL
    shp_ds = gdal.OpenEx(fn_shp, gdal.OF_VECTOR)
    #rasterize
    opts = gdal.RasterizeOptions(burnValues=[1], bands=[1], allTouched=all_touched)
    gdal.Rasterize(raster_ds, shp_ds, options=opts)

    #get rasterized array
    rast_array = raster_ds.GetRasterBand(1).ReadAsArray()
    rast_array = rast_array.astype(float)
    rast_array[rast_array != 1] = np.nan

    #close shp
    shp_ds = None
    raster_ds=None

    return ~np.isnan(rast_array)

def rasterize_feat_shp_ds(shp_ds,layer_name,feat_name,feat_id,raster_ds,all_touched=False,exclude=False):

    if not exclude:
        str_eq="='"
    else:
        str_eq="!='"

    sql_stat='SELECT * FROM '+layer_name+' WHERE '+feat_name+str_eq+feat_id+"'"

    opts=gdal.RasterizeOptions(burnValues=[1],bands=[1],SQLStatement=sql_stat,allTouched=all_touched)
    gdal.Rasterize(raster_ds,shp_ds,options=opts)


def create_mem_raster_on_ref(raster_in):

    ds = gdal.Open(raster_in, gdalconst.GA_ReadOnly)
    match_geotrans = ds.GetGeoTransform()
    match_proj = ds.GetProjection()
    wide = ds.RasterXSize
    high = ds.RasterYSize

    ds= None

    # Create output file
    # array=np.zeros([wide,high])
    ds_out = gdal.GetDriverByName('MEM').Create('', wide, high, 1, gdalconst.GDT_Float32)
    ds_out.SetGeoTransform(match_geotrans)
    ds_out.SetProjection(match_proj)
    band = ds_out.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    band.Fill(-9999, 0)
    # band.WriteArray(array)
    # band.FlushCache()

    return ds_out

def crop_raster_nodata(raster_in,raster_out):

    tmp_dir=create_tmp_dir_for_outfile(raster_out)

    tmp_raster=os.path.join(tmp_dir,'tmp'+os.path.splitext(os.path.basename(raster_in))[0])
    tmp_nodataextent=os.path.join(tmp_dir,'tmp.shp')

    os.system('gdal_calc.py -A '+raster_in+' --outfile '+tmp_raster+' --NoDataValue 0 --calc="1*(A>0)"')

    os.system('gdal_polygonize.py -8 -b 1 -f "ESRI Shapefile" '+tmp_raster+' '+tmp_nodataextent)

    clip_rast_to_shp_sql(raster_in,raster_out,tmp_nodataextent,None)

    remove_tmp_dir_for_outfile(raster_out)

def mask_bitarray(array_in,bit_conf,bool_relate=None):

    mask=np.zeros(np.shape(array_in),dtype=bool)
    mask=operator.or_(mask,array_in==int(bit_conf[0],2))

    if len(bit_conf)>1:
        for i in range(len(bit_conf)-1):
            mask=bool_relate[i](mask,array_in==int(bit_conf[i+1],2))

    return mask

def filter_nanarray(array_in,oper,oper_args,bool_relate=None):

    #short description:
    #oper passes a list of operators based on the operator library. Example: ">" is "operator.gt" (greater than), "=" is "operator.eq", etc..
    #oper_args passes a list for second argument "b" used for each oper[i](array_in,b)
    #bool_relate is the boolean relation between each operation, works the following: bool_relate[1](bool_relate[0](OPER1,OPER2),OPER3) where bool_relate[0] can for instance be operator.and_ or operator.or_...
    #NOTE: if oper and oper_args list N operations, bool_relate takes a list of N-1 boolean operations

    array_out=array_in.astype(float)
    filt=np.zeros(np.shape(array_out),dtype=bool)
    filt=operator.or_(filt,oper[0](array_out,oper_args[0]))

    #if more than one operator
    if len(oper)>1:
        for i in range(len(oper)-1):
            filt=bool_relate[i](filt,oper[i+1](array_out,oper_args[i+1]))

    array_out[filt]=np.NaN

    return array_out, filt

def read_nanarray(raster_in,nodata=None):

    ds=gdal.Open(raster_in,gdalconst.GA_ReadOnly)
    band=ds.GetRasterBand(1)
    if nodata is None:
        nodata=band.GetNoDataValue()
    array=band.ReadAsArray()

    array=array.astype('float32')
    array[array==nodata]=np.NaN

    return array

def update_nanarray(raster_out,array_out,nodata=None):

    ds=gdal.Open(raster_out,gdalconst.GA_Update)
    band=ds.GetRasterBand(1)
    if nodata is None:
        nodata=band.GetNoDataValue()

        array_out[np.isnan(array_out)]=nodata

    band.WriteArray(array_out)

def write_nanarray(raster_out,raster_ref,array_out,nodata=None):

    shutil.copy(raster_ref,raster_out)
    update_nanarray(raster_out,array_out,nodata)

def read_proj_rast(raster_in):

    ds=gdal.Open(raster_in,gdalconst.GA_ReadOnly)
    proj=ds.GetProjection()
    ds=None

    return proj

def pixel_size(fn_raster_in):

    ds = gdal.Open(fn_raster_in, gdal.GA_ReadOnly)
    x0_ref, dx_ref, dxdy_ref, y0_ref, dydx_ref, dy_ref = ds.GetGeoTransform()
    res = np.sqrt(np.abs(dx_ref * dy_ref))

    return res

def extent_rast(raster_in):

    ds = gdal.Open(raster_in, gdalconst.GA_ReadOnly)
    x0_ref, dx_ref, dxdy_ref, y0_ref, dydx_ref, dy_ref = ds.GetGeoTransform()
    proj_wkt = ds.GetProjection()
    col_tot = ds.RasterXSize
    lin_tot = ds.RasterYSize
    x1_ref = x0_ref + col_tot * dx_ref
    y1_ref = y0_ref + lin_tot * dy_ref
    ds = None

    #extent format: Xmin, Ymin, Xmax, Ymax
    xmin=min(x0_ref,x1_ref)
    ymin=min(y0_ref,y1_ref)
    xmax=max(x0_ref,x1_ref)
    ymax=max(y0_ref,y1_ref)

    extent = [xmin, ymin, xmax, ymax]

    return extent, proj_wkt


def translate(raster_in,raster_out,format_out='GTiff',tgt_res=None,interp_method='bilinear'):

    ds=gdal.Open(raster_in)
    opts=gdal.TranslateOptions(format=format_out,xRes=tgt_res[0],yRes=tgt_res[1],resampleAlg=interp_method)
    gdal.Translate(raster_out,ds,options=opts)
    ds=None

def warp_defaultUTM(raster_in,raster_out,format_out='Gtiff',src_EPSG=None,tgt_EPSG=None,tgt_res=None,nodata_in=-9999,nodata_out=-9999,interp_method='bilinear'):

    #default is automatic source proj to UTM of extent centroid, same resolution, with -9999 nodata value and bilinear interpolation
    if tgt_EPSG is None:

        extent, proj_wkt = extent_rast(raster_in)

        poly = poly_from_extent(extent)

        transform = coord_trans_wkt_or_EPSG(True,proj_wkt,False,4326)

        poly.Transform(transform)

        center_lon, center_lat = get_poly_centroid(poly)

        epsg, utm_zone = latlon_to_UTM(center_lat, center_lon)
        print('Projecting tile ' + raster_in + ' to UTM zone: ' + utm_zone + ' corresponding to ESPG: ' + epsg)
        dst_SRS = 'EPSG:'+epsg
    else:
        epsg = str(tgt_EPSG)
        dst_SRS = 'EPSG:'+epsg
        print('Projecting tile ' + raster_in + ' to ESPG: ' + str(tgt_EPSG))

    if tgt_res is None:
        xRes = None
        yRes = None
    else:
        xRes = tgt_res[0]
        yRes = tgt_res[1]

    if src_EPSG is None:
        src_SRS = None
    else:
        src_SRS = 'EPSG:'+str(src_EPSG)

    # warp, see options at: https://gdal.org/python/osgeo.gdal-module.html#Warp
    opts=gdal.WarpOptions(format=format_out,xRes=xRes,yRes=yRes,srcSRS=src_SRS,dstSRS=dst_SRS,srcNodata=nodata_in,dstNodata=nodata_out,resampleAlg=interp_method)
    ds_out = gdal.Warp(raster_out,raster_in,options=opts)

    if format_out == 'MEM':
        # if in-memory datasource, return object
        return ds_out
    else:
        # otherwise, close datasource
        ds_out = None

    #using gdalwarp command line
    # # define target resolution
    # if tgt_res is None:
    #     str_res = ''
    # else:
    #     xres = tgt_res[0]
    # yres = tgt_res[1]
    # str_res = ' -tr ' + str(xres) + ' ' + str(yres)
    # print 'Projecting to resolution: ' + str(xres) + ' ' + str(yres)
    #
    # if srs_EPSG is None:
    #     str_srs = ''
    # else:
    #     str_srs = ' -s_srs EPSG:' + str(srs_EPSG)

    # os.system('gdalwarp -overwrite -r ' + interp_method + str_srs+ ' -t_srs EPSG:' + epsg + str_res+' -dstnodata '+str(nodata)+' -of GTiff ' +raster_in+' '+tif_out)

def reproject_on_ref(ds_in,ds_ref,rast_out,format_out='GTiff',nodata_out=-9999,interp_method='bilinear'):

    #read reference raster georeferencing
    match_proj = ds_ref.GetProjection()
    wide = ds_ref.RasterXSize
    high = ds_ref.RasterYSize
    match_geotrans=ds_ref.GetGeoTransform()

    #read source raster projection
    src_proj = ds_in.GetProjection()

    #create tif output file
    ds_out = gdal.GetDriverByName(format_out).Create(rast_out, wide, high, 1, gdalconst.GDT_Float32)
    ds_out.SetGeoTransform(match_geotrans)
    ds_out.SetProjection(match_proj)
    band = ds_out.GetRasterBand(1)
    band.SetNoDataValue(nodata_out)
    band.Fill(nodata_out, 0)
    band.FlushCache()

    #select interpolation method
    if interp_method == 'bilinear':
        method=gdalconst.GRA_Bilinear
    elif interp_method == 'cubic':
        method=gdalconst.GRA_Cubic
    elif interp_method == 'near':
        method=gdalconst.GRA_NearestNeighbour
    elif interp_method == 'cubicspline':
        method=gdalconst.GRA_CubicSpline
    else:
        sys.exit('Interpolation method not recognized.')

    #reproject
    gdal.ReprojectImage(ds_in, ds_out, src_proj, match_proj, method)

    if format_out == 'MEM':
        #if in-memory datasource, return object
        return ds_out
    else:
        #otherwise, close datasource
        ds_out=None

def reproject_on_ref_fn_tif(raster_in,raster_ref,tif_out):

    ds_in=gdal.Open(raster_in,gdalconst.GA_ReadOnly)
    ds_ref=gdal.Open(raster_ref,gdalconst.GA_ReadOnly)

    reproject_on_ref(ds_in,ds_ref,tif_out)

    ds_in=None
    ds_out=None


def merge_rast_list(list_raster,raster_out,tgt_EPSG=None,nodata_in=None,nodata_out=None):

    #reproject to similar proj if specified
    if tgt_EPSG is not None:

        tmp_dir=create_tmp_dir_for_outfile(raster_out)

        list_raster_proj=[]
        for rast in list_raster:
            print(rast)
            ind = list_raster.index(rast)
            rast_proj=tmp_dir+os.path.splitext(os.path.basename(rast))[0]+'_'+str(ind)+'_proj.tif'
            warp_defaultUTM(rast,rast_proj,tgt_EPSG=tgt_EPSG,nodata_out=nodata_out)
            list_raster_proj.append(rast_proj)

        final_list_raster=list_raster_proj
    else:
        final_list_raster=list_raster

    if nodata_in is None:
        str_nodata_in=''
    else:
        str_nodata_in=' -n '+str(nodata_in)

    if nodata_out is None:
        str_nodata_out=''
    else:
        str_nodata_out=' -a_nodata '+str(nodata_out)

    if os.path.exists(raster_out):
        os.remove(raster_out)

    os.system('gdal_merge.py -o '+raster_out+str_nodata_in+str_nodata_out+' '+" ".join(final_list_raster))

    remove_tmp_dir_for_outfile(raster_out)


def clip_rast_to_shp_sql(raster_in,raster_out,shp_in,sql_req):

    #reproject shp to similar proj if needed
    rast_proj = read_proj_rast(raster_in)
    shp_proj= read_proj_shp(shp_in)

    comp_proj, rast_epsg, _ = compare_proj_wkt_or_EPSG(True,rast_proj,True,shp_proj)

    if not comp_proj:
        tmp_dir=create_tmp_dir_for_outfile(raster_out)

        tmp_shp=tmp_dir+os.path.splitext(os.path.basename(shp_in))[0]+'_tmp.shp'
        os.system('ogr2ogr -t_srs '+rast_epsg+' '+tmp_shp+' '+shp_in)
        final_shp=tmp_shp
    else:
        final_shp=shp_in

    #example of sql request
    #sql_req = "SELECT * FROM your_shapefile"

    gdal.Warp(raster_out,raster_in,cutlineDSName=final_shp,cutlineSQL=sql_req,cropToCutline=True)

    remove_tmp_dir_for_outfile(raster_out)


def clip_rast_to_extent(raster_in,raster_out,extent,extent_EPSG):

    # rast_proj = read_proj_rast(raster_in)
    # comp_proj, _, _ = compare_proj_wkt_or_EPSG(True,rast_proj,False,extent_EPSG)
    #
    # if not comp_proj:
    #     trans=coord_trans_wkt_or_EPSG(False,extent_EPSG,True,rast_proj)
    #
    #     poly=poly_from_extent(extent)
    #
    #     poly.Transform(trans)
    #
    #     final_extent=extent_from_poly(poly)
    # else:
    #     final_extent=extent

    opts=gdal.WarpOptions(outputBounds=(extent[0], extent[1], extent[2], extent[3]),outputBoundsSRS=extent_EPSG)
    gdal.Warp(raster_out,raster_in,options=opts)


#David Shean code
def inters_raster(fn_raster1,fn_raster2):

    def get_ds_srs(ds):
        """Get srs object for GDAL Datset
        """
        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromWkt(ds.GetProjectionRef())
        return ds_srs

    def applyGeoTransform(inX, inY, geoTransform):
        inX = np.asarray(inX)
        inY = np.asarray(inY)
        outX = geoTransform[0] + inX * geoTransform[1] + inY * geoTransform[2]
        outY = geoTransform[3] + inX * geoTransform[4] + inY * geoTransform[5]
        return outX, outY

    def pixelToMap(pX, pY, geoTransform):
        """Convert pixel coordinates to map coordinates based on geotransform

        Accepts float or NumPy arrays
        GDAL model used here - upper left corner of upper left pixel for mX, mY (and in GeoTransform)
        """
        pX = np.asarray(pX, dtype=float)
        pY = np.asarray(pY, dtype=float)
        pX += 0.5
        pY += 0.5
        mX, mY = applyGeoTransform(pX, pY, geoTransform)
        return mX, mY


    def geom_transform(geom, t_srs):
        """Transform a geometry in place
        """
        s_srs = geom.GetSpatialReference()
        if not s_srs.IsSame(t_srs):
            ct = osr.CoordinateTransformation(s_srs, t_srs)
            geom.Transform(ct)
            geom.AssignSpatialReference(t_srs)

    def ds_geom(ds, t_srs=None):
        """Return dataset bbox envelope as geom
        """
        gt = ds.GetGeoTransform()
        ds_srs = get_ds_srs(ds)
        if t_srs is None:
            t_srs = ds_srs
        ns = ds.RasterXSize
        nl = ds.RasterYSize
        x = np.array([0, ns, ns, 0, 0], dtype=float)
        y = np.array([0, 0, nl, nl, 0], dtype=float)
        #Note: pixelToMap adds 0.5 to input coords, need to account for this here
        x -= 0.5
        y -= 0.5
        mx, my = pixelToMap(x, y, gt)
        geom_wkt = 'POLYGON(({0}))'.format(', '.join(['{0} {1}'.format(*a) for a in zip(mx,my)]))
        geom = ogr.CreateGeometryFromWkt(geom_wkt)
        geom.AssignSpatialReference(ds_srs)
        if not ds_srs.IsSame(t_srs):
            geom_transform(geom, t_srs)
        return geom

    def geom_intersection(geom_list, **kwargs):
        convex = False
        intsect = geom_list[0]
        valid = False
        for geom in geom_list[1:]:
            if intsect.Intersects(geom):
                valid = True
                intsect = intsect.Intersection(geom)
        if convex:
            intsect = intsect.ConvexHull()
        if not valid:
            intsect = None
        return intsect

    def ds_geom_intersection(ds_list, **kwargs):
        ref_srs = get_ds_srs(ds_list[0])
        if 't_srs' in kwargs:
            if kwargs['t_srs'] is not None:
                if not ref_srs.IsSame(kwargs['t_srs']):
                    ref_srs = kwargs['t_srs']
        geom_list = []
        for ds in ds_list:
            geom_list.append(ds_geom(ds, t_srs=ref_srs))
        intsect = geom_intersection(geom_list)
        return intsect

    def ds_geom_intersection_extent(ds_list, **kwargs):
        intsect = ds_geom_intersection(ds_list, **kwargs)
        if intsect is not None:
            # Envelope is ul_x, ur_x, lr_y, lr_x
            # Define new geom class with better Envelope options?
            env = intsect.GetEnvelope()
            intsect = [env[0], env[2], env[1], env[3]]
        return intsect

    ds1 = gdal.Open(fn_raster1)
    ds2 = gdal.Open(fn_raster2)

    tag_inters= ds_geom_intersection_extent([ds1,ds2])

    ds1 = None
    ds2 = None

    return tag_inters

