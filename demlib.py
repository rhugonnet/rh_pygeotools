# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

DEM LIBRARY:
Library of Python functions for DEM manipulation: mostly SAGA/ASP wrappers in Python / also some adapted functions from Bob
"""
from __future__ import print_function
import os, sys
import shutil
from shlib import create_tmp_dir_for_outfile, remove_tmp_dir_for_outfile
from osgeo import gdal, osr, ogr, gdalconst

#don't display plots and warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt

#for coreg
import errno
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from pybob.GeoImg import GeoImg
from pybob.ICESat import ICESat
from pybob.coreg_tools import false_hillshade, RMSE, get_geoimg, final_histogram, coreg_fitting, create_stable_mask, get_aspect, get_slope

def dem_slope_mem(fn_dem_in):

    slope = gdal.DEMProcessing('', fn_dem_in, 'slope', format='MEM')

    return slope.GetRasterBand(1).ReadAsArray()

def dem_aspect(fn_dem_in,fn_slope_out):
    pass

def dem_contour_fl(fn_dem_in,fn_shp_out,fl):

    src_ds = gdal.Open(fn_dem_in, gdalconst.GA_ReadOnly)
    src_band = src_ds.GetRasterBand(1)

    mask_band = src_band.GetMaskBand()

    srs = osr.SpatialReference()
    srs.ImportFromWkt(src_ds.GetProjectionRef())

    # remove the .shp if specified by user
    # if fn_shp_out[-4:] == '.shp':
    #     fn_shp_out = fn_shp_out[:-4]

    drv = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(fn_shp_out):
        drv.DeleteDataSource(fn_shp_out)
    dst_ds = drv.CreateDataSource(fn_shp_out)
    dst_layer = dst_ds.CreateLayer('contour', srs=srs)
    field_defn = ogr.FieldDefn('ID', ogr.OFTInteger)
    dst_layer.CreateField(field_defn)
    field_defn = ogr.FieldDefn('elev', ogr.OFTReal)
    dst_layer.CreateField(field_defn)

    #first parameter is interval, third is fixed elevation list to contour
    gdal.ContourGenerate(src_band, 0, 0, fl, 0, 0, dst_layer, 0, 1)
    dst_ds.Destroy()
    src_ds = None


def dem_difference(dem1, dem2, out_dh, glaciermask=None, landmask=None):
    outdir = create_tmp_dir_for_outfile(out_dh)
    master, slave = dem_coregistration_custom(dem1, dem2, glaciermask, landmask, outdir)[0:2]
    master.unmask()
    slave.unmask()

    dh = master.copy(new_raster=(master.img - slave.img))

    dh.write(out_dh)


def dem_mosaic(list_dem_in,dem_out,erode_length=0,hole_fill_length=0,blending_method='mean'):

    # requires asp>2.4
    tmp_mosaic=os.path.dirname(dem_out)+os.path.sep+'tmp'
    os.system('dem_mosaic --erode-length '+ str(erode_length) + ' --hole-fill-length ' + str(hole_fill_length)
              + ' --'+blending_method + ' ' +" ".join(list_dem_in)+' -o '+tmp_mosaic)
    shutil.move(os.path.dirname(dem_out)+os.path.sep + 'tmp-tile-0-'+blending_method+'.tif', dem_out)

def ellipsoid_to_geoid(dem_in,dem_out,geoid='EGM96'):

    #requires asp>2.4
    tmp_geoid=os.path.dirname(dem_out)+os.path.sep+'tmp'
    os.system('dem_geoid '+dem_in+ ' --geoid '+geoid+ ' -o '+tmp_geoid)
    shutil.move(os.path.dirname(dem_out)+os.path.sep+'tmp-adj.tif',dem_out)

def geoid_to_ellipsoid(dem_in,dem_out,geoid='EGM96'):

    #requires asp>2.4
    tmp_ellipsoid=os.path.dirname(dem_out)+os.path.sep+'tmp'
    os.system('dem_geoid '+dem_in+' --geoid '+geoid+' -o '+tmp_ellipsoid+' --reverse-adjustment')
    shutil.move(os.path.dirname(dem_out)+os.path.sep+'tmp-adj.tif',dem_out)
    
def saga_aspect_slope_curvature(dem_in,topo_out,method_nb=5):

    """
    :param dem_in: dem input
    :param topo_out: raster out (3 bands: aspect, slope, curvature)
    :param method_nb: algorithm used, see function description
    :return:

    requirement: SAGA 2.X with X>3

    #ref:http://www.saga-gis.org/saga_tool_doc/2.2.3/ta_morphometry_0.html

    aspect, slope, curvature methods
    methods number
    [0] maximum slope (Travis et al. 1975)
    [1] maximum triangle slope (Tarboton 1997)
    [2] least squares fitted plane (Horn 1981, Costa-Cabral & Burgess 1996)
    [3] 6 parameter 2nd order polynom (Evans 1979)
    [4] 6 parameter 2nd order polynom (Heerdegen & Beran 1982)
    [5] 6 parameter 2nd order polynom (Bauer, Rohdenburg, Bork 1985)
    [6] 9 parameter 2nd order polynom (Zevenbergen & Thorne 1987)
    [7] 10 parameter 3rd order polynom (Haralick 1983)

    unit slope 0=rad, 1=deg, 2=percent
    unit aspect 0=rad, 1=deg
    """

    tmp_dir=create_tmp_dir_for_outfile(topo_out)
    
    saga_elev = tmp_dir + 'elev_temp.sgrd'
    saga_slope = tmp_dir + 'slope_temp.sgrd'
    saga_aspect = tmp_dir + 'aspect_temp.sgrd'
    saga_max_curv = tmp_dir + 'max_curv_temp.sgrd'
    tif_slope = tmp_dir + 'slope_temp.tif'
    tif_aspect = tmp_dir + 'aspect_temp.tif'
    tif_max_curv = tmp_dir + 'max_curv_temp.tif'
    output_vrt = tmp_dir + 'stack.vrt'

    os.system('saga_cmd io_gdal 0 -GRIDS ' + saga_elev + ' -FILES ' + dem_in)
    os.system('saga_cmd ta_morphometry 0 -ELEVATION ' + saga_elev + ' -SLOPE ' + saga_slope + ' -ASPECT ' + saga_aspect + ' -C_MAXI ' + saga_max_curv + ' -METHOD ' +str(method_nb)+' -UNIT_SLOPE 1 -UNIT_ASPECT 1')
    os.system('saga_cmd io_gdal 2 -GRIDS ' + saga_slope + ' -FILE ' + tif_slope)
    os.system('saga_cmd io_gdal 2 -GRIDS ' + saga_aspect + ' -FILE ' + tif_aspect)
    os.system('saga_cmd io_gdal 2 -GRIDS ' + saga_max_curv + ' -FILE ' + tif_max_curv)

    os.system('gdalbuildvrt -separate -overwrite ' + output_vrt + ' ' + tif_slope + ' ' + tif_aspect + ' ' + tif_max_curv)
    os.system('gdal_translate ' + output_vrt + ' ' + topo_out)

    remove_tmp_dir_for_outfile(topo_out)


def dem_coregistration_custom(masterDEM, slaveDEM, glaciermask=None, landmask=None, outdir='.', pts=False, full_ext=False,magnlimit=1.):
    #This coreg code is from Robert McNaab and Chris Nuth: pybob
    """
    Iteratively co-register elevation data, based on routines described in Nuth and Kaeaeb, 2011.

    Parameters
    ----------
    masterDEM : string or GeoImg
        Path to filename or GeoImg dataset representing "master" DEM.
    slaveDEM : string or GeoImg
        Path to filename or GeoImg dataset representing "slave" DEM.
    glaciermask : string, optional
        Path to shapefile representing points to exclude from co-registration
        consideration (i.e., glaciers).
    landmask : string, optional
        Path to shapefile representing points to include in co-registration
        consideration (i.e., stable ground/land).
    outdir : string, optional
        Location to save co-registration outputs.
    pts : bool, optional
        If True, program assumes that masterDEM represents point data (i.e., ICESat),
        as opposed to raster data. Slope/aspect are then calculated from slaveDEM.
        masterDEM should be a string representing an HDF5 file continaing ICESat data.
    full_ext : bool, optional
        If True, program writes full extents of input DEMs. If False, program writes
        input DEMs cropped to their common extent. Default is False.
    """

    def preprocess(stable_mask, slope, aspect, master, slave):

        if isinstance(master, GeoImg):
            stan = np.tan(np.radians(slope)).astype(np.float32)
            dH = master.copy(new_raster=(master.img - slave.img))
            dH.img[stable_mask] = np.nan
            master_mask = isinstance(master.img, np.ma.masked_array)
            slave_mask = isinstance(slave.img, np.ma.masked_array)

            if master_mask and slave_mask:
                dH.mask(np.logical_or(master.img.mask, slave.img.mask))
            elif master_mask:
                dH.mask(master.img.mask)
            elif slave_mask:
                dH.mask(slave.img.mask)

            if dH.isfloat:
                dH.img[stable_mask] = np.nan

            #adding a 5NMAD filtering for robustness at various resolutions
            myfirstkeep = ((np.absolute(dH.img) < 200.0) & np.isfinite(dH.img) & (aspect>0))
            nmad = 1.4826 * np.median(np.abs(dH.img[myfirstkeep]) - np.median(dH.img[myfirstkeep]))

            dHtan = dH.img / stan
            #here too
            mykeep = ((np.absolute(dH.img) < 200.0) & np.isfinite(dH.img) &
                      (slope > 7.0) & (dH.img != 0.0) & (aspect >= 0) & (np.absolute(dH.img - np.median(dH.img[myfirstkeep])) < 5*nmad))
            dH.img[np.invert(mykeep)] = np.nan
            xdata = aspect[mykeep]
            ydata = dHtan[mykeep]
            sdata = stan[mykeep]

        elif isinstance(master, ICESat):
            slave_pts = slave.raster_points(master.xy)
            dH = master.elev - slave_pts

            slope_pts = slope.raster_points(master.xy)
            stan = np.tan(np.radians(slope_pts))

            aspect_pts = aspect.raster_points(master.xy)
            smask = stable_mask.raster_points(master.xy) > 0

            dH[smask] = np.nan

            dHtan = dH / stan

            mykeep = ((np.absolute(dH) < 200.0) & np.isfinite(dH) &
                      (slope_pts > 3.0) & (dH != 0.0) & (aspect_pts >= 0))

            dH[np.invert(mykeep)] = np.nan
            xdata = aspect_pts[mykeep]
            ydata = dHtan[mykeep]
            sdata = stan[mykeep]

        return dH, xdata, ydata, sdata


    # if the output directory does not exist, create it.
    outdir = os.path.abspath(outdir)
    try:
        os.makedirs(outdir)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else:
            raise
    # make a file to save the coregistration parameters to.
    paramf = open(outdir + os.path.sep + 'coreg_params.txt', 'w')
    # create the output pdf
    pp = PdfPages(outdir + os.path.sep + 'CoRegistration_Results.pdf')

    if full_ext:
        print('Writing full extents of output DEMs.')
    else:
        print('Writing DEMs cropped to common extent.')

    if type(masterDEM) is str:
        mfilename = os.path.basename(masterDEM)
        mfiledir = os.path.dirname(masterDEM)
    else:
        mfilename = masterDEM.filename
        mfiledir = masterDEM.in_dir_path

    if type(slaveDEM) is str:
        sfilename = os.path.basename(slaveDEM)
    else:
        sfilename = slaveDEM.filename

    slaveDEM = get_geoimg(slaveDEM)
    # if we're dealing with ICESat/pt data, change how we load masterDEM data
    if pts:
        masterDEM = ICESat(masterDEM)
        masterDEM.project('epsg:{}'.format(slaveDEM.epsg))
        mybounds = [slaveDEM.xmin, slaveDEM.xmax, slaveDEM.ymin, slaveDEM.ymax]
        masterDEM.clip(mybounds)
        masterDEM.clean()
        slope_geo = get_slope(slaveDEM)
        aspect_geo = get_aspect(slaveDEM)

        slope_geo.write('tmp_slope.tif', out_folder=outdir)
        aspect_geo.write('tmp_aspect.tif', out_folder=outdir)

        smask = create_stable_mask(slaveDEM, glaciermask, landmask)
        slaveDEM.mask(smask)
        stable_mask = slaveDEM.copy(new_raster=smask)  # make the mask a geoimg
    else:
        orig_masterDEM = get_geoimg(masterDEM)

        masterDEM = orig_masterDEM.reproject(slaveDEM)  # need to resample masterDEM to cell size of slave.
        # masterDEM.img[masterDEM.img<1]=np.nan
        stable_mask = create_stable_mask(masterDEM, glaciermask, landmask)

        slope_geo = get_slope(masterDEM)
        aspect_geo = get_aspect(masterDEM)
        slope_geo.write('tmp_slope.tif', out_folder=outdir)
        aspect_geo.write('tmp_aspect.tif', out_folder=outdir)
        masterDEM.mask(stable_mask)

    slope = slope_geo.img
    aspect = aspect_geo.img

    mythresh = np.float64(200)  # float64 really necessary?
    mystd = np.float64(200)
    mycount = 0
    tot_dx = np.float64(0)
    tot_dy = np.float64(0)
    tot_dz = np.float64(0)
    magnthresh = 200
    # magnlimit = 1
    mytitle = 'DEM difference: pre-coregistration'
    if pts:
        this_slave = slaveDEM
        this_slave.mask(stable_mask.img)
    else:
        this_slave = slaveDEM.reproject(masterDEM)
        this_slave.mask(stable_mask)

    while mythresh > 2 and magnthresh > magnlimit:
        if mycount != 0:
            # slaves.append(slaves[-1].reproject(masterDEM))
            # slaves[-1].mask(stable_mask)
            mytitle = "DEM difference: After Iteration {}".format(mycount)
        mycount += 1
        print("Running iteration #{}".format(mycount))
        print("Running iteration #{}".format(mycount), file=paramf)

        # if we don't have two DEMs, showing the false hillshade doesn't work.
        if not pts:
            dH, xdata, ydata, sdata = preprocess(stable_mask, slope, aspect, masterDEM, this_slave)
            false_hillshade(dH, mytitle, pp)
            dH_img = dH.img
        else:
            dH, xdata, ydata, sdata = preprocess(stable_mask, slope_geo, aspect_geo, masterDEM, this_slave)
            dH_img = dH

        if mycount == 1:
            dH0 = dH_img

        # calculate threshold, standard deviation of dH
        # mythresh = 100 * (mystd-np.nanstd(dH_img))/mystd
        # mystd = np.nanstd(dH_img)
        # USE RMSE instead ( this is to make su that there is improvement in the spread)
        mythresh = 100 * (mystd - RMSE(dH_img)) / mystd
        mystd = RMSE(dH_img)

        mytitle2 = "Co-registration: Iteration {}".format(mycount)
        dx, dy, dz = coreg_fitting(xdata, ydata, sdata, mytitle2, pp)
        tot_dx += dx
        tot_dy += dy
        tot_dz += dz
        magnthresh = np.sqrt(np.square(dx) + np.square(dy) + np.square(dz))
        print(tot_dx, tot_dy, tot_dz)
        print(tot_dx, tot_dy, tot_dz, file=paramf)
        # print np.nanmean(slaves[-1].img)

        # print slaves[-1].xmin, slaves[-1].ymin

        # shift most recent slave DEM
        this_slave.shift(dx, dy)  # shift in x,y
        # print tot_dx, tot_dy
        # no idea why slaves[-1].img += dz doesn't work, but the below seems to.
        zupdate = np.ma.array(this_slave.img.data + dz, mask=this_slave.img.mask)  # shift in z
        this_slave = this_slave.copy(new_raster=zupdate)
        if pts:
            this_slave.mask(stable_mask.img)
            slope_geo.shift(dx, dy)
            aspect_geo.shift(dx, dy)
            stable_mask.shift(dx, dy)
        else:
            this_slave = this_slave.reproject(masterDEM)
            this_slave.mask(stable_mask)

        print("Percent-improvement threshold and Magnitute threshold:")
        print(mythresh, magnthresh)

        # slaves[-1].display()
        if mythresh > 2 and magnthresh > magnlimit:
            dH = None
            dx = None
            dy = None
            dz = None
            xdata = None
            ydata = None
            sdata = None
        else:
            if not pts:
                dH, xdata, ydata, sdata = preprocess(stable_mask, slope, aspect, masterDEM, this_slave)
                mytitle = "DEM difference: After Iteration {}".format(mycount)
                # adjust final dH
                # myfadj=np.nanmean([np.nanmean(dH.img),np.nanmedian(dH.img)])
                # myfadj=np.nanmedian(dH.img)
                # tot_dz += myfadj
                # dH.img = dH.img-myfadj

                false_hillshade(dH, mytitle, pp)
                dHfinal = dH.img
            else:
                mytitle2 = "Co-registration: FINAL"
                dH, xdata, ydata, sdata = preprocess(stable_mask, slope_geo, aspect_geo, masterDEM, this_slave)
                dx, dy, dz = coreg_fitting(xdata, ydata, sdata, mytitle2, pp)
                dHfinal = dH

    # Create final histograms pre and post coregistration
    # shift = [tot_dx, tot_dy, tot_dz]  # commented because it wasn't actually used.
    final_histogram(dH0, dHfinal, pp)

    # create new raster with dH sample used for co-registration as the band
    # dHSample = dH.copy(new_raster=dHpost_sample)
    # dHSample.write(outdir + os.path.sep + 'dHpost_sample.tif') # have to fill these in!
    # save full dH output
    # dHfinal.write('dHpost.tif', out_folder=outdir)
    # save adjusted slave dem
    if sfilename is not None:
        slaveoutfile = '.'.join(sfilename.split('.')[0:-1]) + '_adj.tif'
    else:
        slaveoutfile = 'slave_adj.tif'

    if pts:
        outslave = slaveDEM.copy()
    else:
        if full_ext:
            outslave = get_geoimg(slaveDEM)
        else:
            outslave = slaveDEM.reproject(masterDEM)

    outslave.shift(tot_dx, tot_dy)
    outslave.img = outslave.img + tot_dz
    outslave.write(slaveoutfile, out_folder=outdir)
    outslave.filename = slaveoutfile

    if pts:
        slope_geo.write('tmp_slope.tif', out_folder=outdir)
        aspect_geo.write('tmp_aspect.tif', out_folder=outdir)

    # Final Check --- for debug
    if not pts:
        dH, xdata, ydata, sdata = preprocess(stable_mask, slope, aspect, masterDEM, outslave)
        false_hillshade(dH, 'FINAL CHECK', pp)

        if mfilename is not None:
            mastoutfile = '.'.join(mfilename.split('.')[0:-1]) + '_adj.tif'
        else:
            mastoutfile = 'master_adj.tif'

        if full_ext:
            masterDEM = orig_masterDEM
        masterDEM.write(mastoutfile, out_folder=outdir)

    pp.close()
    print("Fin.")
    print("Fin.", file=paramf)
    paramf.close()
    plt.close('all')

    out_offs = [tot_dx, tot_dy, tot_dz]

    return masterDEM, outslave, out_offs, mystd

def valid_points_coreg(masterDEM, slaveDEM, glaciermask=None, landmask=None, pts=False):

    if type(masterDEM) is str:
        mfilename = os.path.basename(masterDEM)
        mfiledir = os.path.dirname(masterDEM)
    else:
        mfilename = masterDEM.filename
        mfiledir = masterDEM.in_dir_path

    if type(slaveDEM) is str:
        sfilename = os.path.basename(slaveDEM)
    else:
        sfilename = slaveDEM.filename

    slaveDEM = get_geoimg(slaveDEM)
    # if we're dealing with ICESat/pt data, change how we load masterDEM data
    if pts:
        masterDEM = ICESat(masterDEM)
        masterDEM.project('epsg:{}'.format(slaveDEM.epsg))
        mybounds = [slaveDEM.xmin, slaveDEM.xmax, slaveDEM.ymin, slaveDEM.ymax]
        masterDEM.clip(mybounds)
        masterDEM.clean()
        slope_geo = get_slope(slaveDEM)
        aspect_geo = get_aspect(slaveDEM)

        smask = create_stable_mask(slaveDEM, glaciermask, landmask)
        slaveDEM.mask(smask)
        stable_mask = slaveDEM.copy(new_raster=smask)  # make the mask a geoimg
    else:
        orig_masterDEM = get_geoimg(masterDEM)

        masterDEM = orig_masterDEM.reproject(slaveDEM)  # need to resample masterDEM to cell size of slave.
        # masterDEM.img[masterDEM.img<1]=np.nan
        stable_mask = create_stable_mask(masterDEM, glaciermask, landmask)

        slope_geo = get_slope(masterDEM)
        aspect_geo = get_aspect(masterDEM)
        masterDEM.mask(stable_mask)

    slope = slope_geo.img
    aspect = aspect_geo.img

    if pts:
        this_slave = slaveDEM
        this_slave.mask(stable_mask.img)
    else:
        this_slave = slaveDEM.reproject(masterDEM)
        this_slave.mask(stable_mask)

    master=masterDEM
    slave=this_slave

    if isinstance(master, GeoImg):
        stan = np.tan(np.radians(slope)).astype(np.float32)
        dH = master.copy(new_raster=(master.img - slave.img))
        dH.img[stable_mask] = np.nan
        master_mask = isinstance(master.img, np.ma.masked_array)
        slave_mask = isinstance(slave.img, np.ma.masked_array)

        if master_mask and slave_mask:
            dH.mask(np.logical_or(master.img.mask, slave.img.mask))
        elif master_mask:
            dH.mask(master.img.mask)
        elif slave_mask:
            dH.mask(slave.img.mask)

        if dH.isfloat:
            dH.img[stable_mask] = np.nan

        # adding a 5NMAD filtering for robustness at various resolutions
        myfirstkeep = ((np.absolute(dH.img) < 200.0) & np.isfinite(dH.img) & (aspect > 0))
        nmad = 1.4826 * np.median(np.abs(dH.img[myfirstkeep]) - np.median(dH.img[myfirstkeep]))

        dHtan = dH.img / stan
        # here too
        mykeep = ((np.absolute(dH.img) < 200.0) & np.isfinite(dH.img) &
                  (slope > 7.0) & (dH.img != 0.0) & (aspect >= 0) & (
                              np.absolute(dH.img - np.median(dH.img[myfirstkeep])) < 5 * nmad))
        dH.img[np.invert(mykeep)] = np.nan

    elif isinstance(master, ICESat):
        slave_pts = slave.raster_points(master.xy)
        dH = master.elev - slave_pts

        slope_pts = slope.raster_points(master.xy)
        stan = np.tan(np.radians(slope_pts))

        aspect_pts = aspect.raster_points(master.xy)
        smask = stable_mask.raster_points(master.xy) > 0

        dH[smask] = np.nan

        dHtan = dH / stan

        mykeep = ((np.absolute(dH) < 200.0) & np.isfinite(dH) &
                  (slope_pts > 3.0) & (dH != 0.0) & (aspect_pts >= 0))

        dH[np.invert(mykeep)] = np.nan
    else:
        sys.exit('Master DEM not recognized.')

    return np.count_nonzero(mykeep)
