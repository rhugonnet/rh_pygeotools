# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

DEM STRIP PRODUCTS LIBRARY:
Library of Python functions to process DEM strip products such as:
- ArcticDEM release 7 strip (2m)
- REMA release 1.1 strip (2m and 8m)
- ASTER L1A via ASP
- ASTER L1A via MMASTERv2

All strips placed within the input directory are processed to a specified output directory

3 different options for output directory structure:
1/ all strips at the root of directory
2/ all strips in subdirectories following input structure
3/ all strips are referenced to a 1x1° lat/lon tile subdirectory based on their centroid


"""
from __future__ import print_function
import os, sys, shutil, csv
from shlib import create_tmp_dir_for_outfile, remove_tmp_dir_for_outfile, extract_file_from_tar_gz
from rastlib import warp_defaultUTM, inters_raster, translate
from demlib import ellipsoid_to_geoid, geoid_to_ellipsoid, dem_coregistration_custom, dem_mosaic, valid_points_coreg


def ArcticDEM_strip_r7(dir_ArcDEM,tile_name,dem_dir_out,mosaic_coreg_segm=True,filter_params=None,format_out='GTiff',tgt_EPSG=3413,tgt_res=None,nodata_out=-9999,interp_method=None,geoid=False,rm_tar=False,downsample=True):

    print('--ArcticDEM release 7 strip processing--')
    # 1/ LOCATE STRIPS

    print('Searching for tile ' + tile_name + ' in folder '+ dir_ArcDEM +'...')
    subtile_dir = os.path.join(dir_ArcDEM, tile_name)
    seg_tar_gz_list = [os.path.join(subtile_dir,tar_file) for tar_file in os.listdir(subtile_dir) if
                        tar_file.endswith('.tar.gz')]
    print('Found '+str(len(seg_tar_gz_list))+' segments in tile folder.')

    # 2/ EXTRACT ALL STRIPS

    tmp_dir = create_tmp_dir_for_outfile(os.path.join(dem_dir_out,'all_strips'))

    list_tmp_dem = [os.path.join(tmp_dir, os.path.splitext(os.path.splitext(os.path.basename(seg_tar_gz))[0])[0] + '_dem.tif') for seg_tar_gz in
                    seg_tar_gz_list]
    for seg_tar_gz in seg_tar_gz_list:
        print('Extracting dem file of segment '+str(seg_tar_gz_list.index(seg_tar_gz)+1)+' out of '+str(len(seg_tar_gz_list)))
        extract_file_from_tar_gz(seg_tar_gz, os.path.splitext(os.path.splitext(os.path.basename(seg_tar_gz))[0])[0] + '_dem.tif',
                                 list_tmp_dem[seg_tar_gz_list.index(seg_tar_gz)])

    if rm_tar:
        #once extracted, remove original data to avoid taking up too much space
        shutil.rmtree(subtile_dir)

    if downsample:
        list_ds=[]
        for tmp_dem in list_tmp_dem:
            tmp_dem_ds = os.path.join(tmp_dir, os.path.splitext(os.path.basename(tmp_dem))[0] + '_ds.tif')
            list_ds.append(tmp_dem_ds)
            translate(tmp_dem, tmp_dem_ds, format_out='GTiff', tgt_res=[8,-8], interp_method='bilinear')
            os.remove(tmp_dem)
        list_tmp_dem=list_ds

    # 2.5/ MOSAIC/COREGISTER SETSM SEGMENTS ITERATIVELY INTO STRIPS

    if mosaic_coreg_segm:

        max_coreg_rmse = 10.
        nb_init_pixel_min = 1000

        if not os.path.exists(os.path.join(dem_dir_out, 'coreg_results')):
            os.mkdir(os.path.join(dem_dir_out, 'coreg_results'))
        if not os.path.exists(os.path.join(dem_dir_out, 'coreg_results', 'good')):
            os.mkdir(os.path.join(dem_dir_out, 'coreg_results', 'good'))
        if not os.path.exists(os.path.join(dem_dir_out, 'coreg_results', 'breaks')):
            os.mkdir(os.path.join(dem_dir_out, 'coreg_results', 'breaks'))

        summary = os.path.join(dem_dir_out, 'mosaic_coreg_summary_'+tile_name+'.csv')

        with open(summary,'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['scene','segment','master','dx','dy','dz','RMSE','mosaic_nb'])

            #create segment lists
            list_seg = list(list_tmp_dem)
            list_seg_id = []
            for seg in list_seg:
                seg_id_split = os.path.basename(seg).split('_')
                seg_id='_'.join(seg_id_split[0:5])
                seg_suf = '_'.join(seg_id_split[6:])
                # seg_nb=int(seg_id_split[5][3:])

                if seg_id not in list_seg_id:
                    list_seg_id.append(seg_id)

            #list all list of scenes acquired the same day
            list_list_by_seg_id=[]
            for seg_id in list_seg_id:
                tmp_list = [seg for seg in list_seg if seg_id in seg]
                tmp_list = sorted(tmp_list, key = lambda x: x.split('/')[-1].split('_')[5][3:].zfill(2))

                list_list_by_seg_id.append(tmp_list)

            # coregister in series
            list_final_strip=[]

            print('List of list of segments:')
            print(list_list_by_seg_id)

            #loop through same day scenes
            for list_by_seg_id in list_list_by_seg_id:

                print(list_by_seg_id)

                #bother to coreg only if there is more than 1 segment
                if len(list_by_seg_id)>1:

                    tmp_dir_coreg = create_tmp_dir_for_outfile(os.path.join(tmp_dir,'coreg'))
                    #master dem is the first segment
                    dem_master = list_by_seg_id[0]
                    k=0

                    #for csv summary: first segment
                    dem_master_id_split = os.path.basename(dem_master).split('_')
                    dem_master_id = '_'.join(dem_master_id_split[0:5])
                    dem_master_seg_nb = int(dem_master_id_split[5][3:])
                    writer.writerow([dem_master_id, dem_master_seg_nb, dem_master_seg_nb, 'first_seg', 'first_seg', 'first_seg', 'first_seg', str(k+1)])

                    for i in range(len(list_by_seg_id)-1):

                        dem_in = list_by_seg_id[i+1]
                        print('Master segment is ' + dem_master)
                        print('Slave segment is '+dem_in)

                        #for csv summary: all segments
                        dem_in_id_split = os.path.basename(dem_in).split('_')
                        dem_in_id = '_'.join(dem_in_id_split[0:5])
                        dem_in_seg_nb = int(dem_in_id_split[5][3:])
                        dem_master_id_split = os.path.basename(dem_master).split('_')
                        dem_master_seg_nb = int(dem_master_id_split[5][3:])
                        k0=k

                        #deriving intersection
                        extent = inters_raster(dem_in,dem_master)

                        nb_init_pixel = valid_points_coreg(dem_master,dem_in)

                        #try coreg only if intersection is non void
                        if extent is not None and nb_init_pixel>nb_init_pixel_min:

                            tmp_dir_seg = create_tmp_dir_for_outfile(dem_in)

                            #do coreg, save results
                            try:
                                out_offs, final_rmse = dem_coregistration_custom(dem_master,dem_in,outdir=tmp_dir_seg)[2:4]
                            except Exception:
                                out_offs = None
                                final_rmse=999

                            if final_rmse < max_coreg_rmse:
                                str_out = 'good'
                            else:
                                str_out = 'breaks'

                            shutil.move(os.path.join(tmp_dir_seg, 'CoRegistration_Results.pdf'),
                                        os.path.join(dem_dir_out, 'coreg_results',str_out,
                                        os.path.splitext(os.path.basename(dem_in))[0] + '_' +
                                        os.path.splitext(os.path.basename(dem_master))[0] + '_coreg_results.pdf'))
                            shutil.move(os.path.join(tmp_dir_seg, 'coreg_params.txt'),
                                        os.path.join(dem_dir_out, 'coreg_results',str_out,
                                        os.path.splitext(os.path.basename(dem_in))[0] + '_' +
                                        os.path.splitext(os.path.basename(dem_master))[0] + '_coreg_params.txt'))

                            print('Final coregistration RMSE is:'+str(final_rmse))

                            #if rmse is good, keep doing mosaics with next segment: mosaic is new master
                            if final_rmse < max_coreg_rmse:
                                print('All good. Doing mosaic and iterating to next segment...')
                                tmp_mosaic_list = [dem_master,os.path.join(tmp_dir_seg,os.path.splitext(os.path.basename(dem_in))[0]+'_adj.tif')]

                                tmp_dem_coreg=os.path.join(tmp_dir_coreg,os.path.basename(list_by_seg_id[i+1]))

                                dem_mosaic(tmp_mosaic_list,tmp_dem_coreg,erode_length=0,hole_fill_length=0,blending_method='mean')
                                dem_master = tmp_dem_coreg

                                #if this is the last segment, save mosaic and add to list of finalized mosaics
                                if i==(len(list_by_seg_id)-2):
                                    k=k+1
                                    final_strip = os.path.join(tmp_dir,
                                                               list_seg_id[list_list_by_seg_id.index(
                                                                   list_by_seg_id)] + '_mosaic' + str(k) + '_' + seg_suf)
                                    shutil.copyfile(dem_master,final_strip)
                                    list_final_strip.append(final_strip)

                            #otherwise, break the iterative mosaic-ing
                            else:
                                #if this is the last segment
                                if i==(len(list_by_seg_id)-2):
                                    #if master is the previous segment alone
                                    if dem_master == list_by_seg_id[i]:
                                        print('RMSE is larger than limit of '+str(max_coreg_rmse)+'m: saving last segments independently...')
                                        list_final_strip.append(dem_in)
                                        list_final_strip.append(dem_master)
                                    #if master is a previous ongoing mosaic
                                    else:
                                        k = k + 1
                                        final_strip = os.path.join(tmp_dir,
                                                                   list_seg_id[list_list_by_seg_id.index(
                                                                       list_by_seg_id)] + '_mosaic' + str(
                                                                       k) + '_' + seg_suf)
                                        shutil.copyfile(dem_master, final_strip)
                                        print('RMSE is larger than limit of ' + str(
                                            max_coreg_rmse) + 'm: saving previous mosaic and last segment...')
                                        list_final_strip.append(dem_in)
                                        list_final_strip.append(final_strip)
                                #if there is still segments to go
                                else:
                                    #save current mosaic
                                    k = k + 1
                                    final_strip = os.path.join(tmp_dir,
                                                               list_seg_id[list_list_by_seg_id.index(
                                                                   list_by_seg_id)] + '_mosaic' + str(
                                                                   k) + '_' + seg_suf)
                                    shutil.copyfile(dem_master, final_strip)
                                    print('RMSE is larger than limit of ' + str(
                                        max_coreg_rmse) + 'm: saving previous mosaic, slave is new master...')
                                    list_final_strip.append(final_strip)
                                    #new master is the segment that did not mosaic
                                    dem_master = dem_in

                            remove_tmp_dir_for_outfile(dem_in)

                        #if no intersection and last segment
                        elif i==(len(list_by_seg_id)-2):
                            out_offs = None
                            final_rmse = None
                            if dem_master == list_by_seg_id[i]:
                                print('No intersection between master and slave: saving last segments independently...')
                                list_final_strip.append(dem_in)
                                list_final_strip.append(dem_master)
                            # if master is a previous ongoing mosaic
                            else:
                                k = k + 1
                                final_strip = os.path.join(tmp_dir,
                                                           list_seg_id[list_list_by_seg_id.index(
                                                               list_by_seg_id)] + '_mosaic' + str(
                                                               k) + '_' + seg_suf)
                                shutil.copyfile(dem_master, final_strip)
                                print('No intersection between master and slave: saving previous mosaic and last segment...')
                                list_final_strip.append(dem_in)
                                list_final_strip.append(final_strip)
                        #if no intersection and still segments to go
                        else:
                            out_offs = None
                            final_rmse = None
                            # save current mosaic
                            k = k + 1
                            final_strip = os.path.join(tmp_dir,
                                                       list_seg_id[list_list_by_seg_id.index(
                                                           list_by_seg_id)] + '_mosaic' + str(
                                                           k) + '_' + seg_suf)
                            shutil.copyfile(dem_master, final_strip)
                            print('No intersection between master and slave: saving previous mosaic, slave is new master...')
                            list_final_strip.append(final_strip)
                            # new master is the segment that did not mosaic
                            dem_master = dem_in

                        #write data in summary
                        if out_offs is None:
                            dx = 'no_inters'
                            dy = 'no_inters'
                            dz = 'no_inters'
                            rmse = 'no_inters'
                        else:
                            dx = str(out_offs[0])
                            dy = str(out_offs[1])
                            dz = str(out_offs[2])
                            rmse = str(final_rmse)

                        writer.writerow([dem_in_id, dem_in_seg_nb, dem_master_seg_nb, dx, dy, dz, rmse, str(k0+1)])


                    remove_tmp_dir_for_outfile(os.path.join(tmp_dir,'coreg'))

                else:
                    list_final_strip.append(list_by_seg_id[0])

                    dem_master_id_split = os.path.basename(list_by_seg_id[0]).split('_')
                    dem_master_id = '_'.join(dem_master_id_split[0:5])
                    dem_master_seg_nb = int(dem_master_id_split[5][3:])
                    writer.writerow(
                        [dem_master_id, dem_master_seg_nb, dem_master_seg_nb, 'only_seg', 'only_seg', 'only_seg',
                         'only_seg', 'only_seg'])


    else:
        list_final_strip = list_tmp_dem


    # 3/ PROCESSING OF FINAL STRIPS

    final_dem = []
    for final_strip in list_final_strip:

        # 3.1/ FILTER STRIPS
        if filter_params is not None:
            sys.exit('No filter pre-defined for this DEM product.')

        # 3.2/ REPROJECT STRIPS

        # raw data is GeoTiff, 3413, 1 arc-sec and -9999 nodata_out
        if format_out == 'GTiff' and tgt_EPSG == 3413 and tgt_res is None and nodata_out is -9999:
            tmp_dem_proj = final_strip
        else:
            tmp_dem_proj = os.path.join(tmp_dir, os.path.splitext(os.path.basename(final_strip))[0] + '_proj.tif')
            warp_defaultUTM(final_strip, tmp_dem_proj, format_out, 3413, tgt_EPSG, tgt_res, nodata_out, interp_method)

        # 3.3/ ELLIPSOID OR GEOID STRIPS

        # raw data is ellipsoid WGS84
        if geoid:
            tmp_dem_geoid = os.path.join(tmp_dir, os.path.splitext(os.path.basename(final_strip))[0] + '_geoid.tif')
            ellipsoid_to_geoid(tmp_dem_proj, tmp_dem_geoid)
        else:
            tmp_dem_geoid = tmp_dem_proj

        final_dem.append(tmp_dem_geoid)

    for dem in final_dem:
        shutil.copyfile(dem,os.path.join(dem_dir_out,os.path.basename(dem)))

    remove_tmp_dir_for_outfile(os.path.join(dem_dir_out,'all_strips'))

    print('Fin.')

def REMA_strip_r1_1(dir_REMA,tile_name,dem_dir_out,mosaic_coreg_segm=True,filter_params=None,format_out='GTiff',tgt_EPSG=3031,tgt_res=None,nodata_out=-9999,interp_method=None,geoid=False,rm_tar=False,downsample=True):

    print('--REMA release 1.1 strip processing--')

    # 1/ LOCATE STRIPS

    print('Searching for tile ' + tile_name + ' in folder ' + dir_REMA + '...')
    subtile_dir = os.path.join(dir_REMA, tile_name)
    seg_tar_gz_list = [os.path.join(subtile_dir, tar_file) for tar_file in os.listdir(subtile_dir) if
                       tar_file.endswith('.tar.gz')]
    print('Found ' + str(len(seg_tar_gz_list)) + ' segments in tile folder.')

    # 2/ EXTRACT ALL STRIPS

    tmp_dir = create_tmp_dir_for_outfile(os.path.join(dem_dir_out, 'all_strips'))

    list_tmp_dem = [
        os.path.join(tmp_dir, os.path.splitext(os.path.splitext(os.path.basename(seg_tar_gz))[0])[0] + '_dem.tif') for
        seg_tar_gz in
        seg_tar_gz_list]
    for seg_tar_gz in seg_tar_gz_list:
        print('Extracting dem file of segment ' + str(seg_tar_gz_list.index(seg_tar_gz) + 1) + ' out of ' + str(
            len(seg_tar_gz_list)))
        extract_file_from_tar_gz(seg_tar_gz,
                                 os.path.splitext(os.path.splitext(os.path.basename(seg_tar_gz))[0])[0] + '_dem.tif',
                                 list_tmp_dem[seg_tar_gz_list.index(seg_tar_gz)])

    if rm_tar:
        # once extracted, remove original data to avoid taking up too much space
        shutil.rmtree(subtile_dir)

    if downsample:
        list_ds = []
        for tmp_dem in list_tmp_dem:
            tmp_dem_ds = os.path.join(tmp_dir, os.path.splitext(os.path.basename(tmp_dem))[0] + '_ds.tif')
            list_ds.append(tmp_dem_ds)
            translate(tmp_dem, tmp_dem_ds, format_out='GTiff', tgt_res=[8, -8], interp_method='bilinear')
            os.remove(tmp_dem)
        list_tmp_dem = list_ds

    # 2.5/ MOSAIC/COREGISTER SETSM SEGMENTS ITERATIVELY INTO STRIPS

    if mosaic_coreg_segm:

        max_coreg_rmse = 10.
        nb_init_pixel_min = 1000

        if not os.path.exists(os.path.join(dem_dir_out, 'coreg_results')):
            os.mkdir(os.path.join(dem_dir_out, 'coreg_results'))
        if not os.path.exists(os.path.join(dem_dir_out, 'coreg_results', 'good')):
            os.mkdir(os.path.join(dem_dir_out, 'coreg_results', 'good'))
        if not os.path.exists(os.path.join(dem_dir_out, 'coreg_results', 'breaks')):
            os.mkdir(os.path.join(dem_dir_out, 'coreg_results', 'breaks'))

        summary = os.path.join(dem_dir_out, 'mosaic_coreg_summary_' + tile_name + '.csv')

        with open(summary, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(['scene', 'segment', 'master', 'dx', 'dy', 'dz', 'RMSE', 'mosaic_nb'])

            # create segment lists
            list_seg = list(list_tmp_dem)
            list_seg_id = []
            for seg in list_seg:
                seg_id_split = os.path.basename(seg).split('_')
                seg_id = '_'.join(seg_id_split[0:5])
                seg_suf = '_'.join(seg_id_split[6:])
                # seg_nb=int(seg_id_split[5][3:])

                if seg_id not in list_seg_id:
                    list_seg_id.append(seg_id)

            # list all list of scenes acquired the same day
            list_list_by_seg_id = []
            for seg_id in list_seg_id:
                tmp_list = [seg for seg in list_seg if seg_id in seg]
                tmp_list = sorted(tmp_list, key=lambda x: x.split('/')[-1].split('_')[5][3:].zfill(2))

                list_list_by_seg_id.append(tmp_list)

            # coregister in series
            list_final_strip = []

            print('List of list of segments:')
            print(list_list_by_seg_id)

            # loop through same day scenes
            for list_by_seg_id in list_list_by_seg_id:

                print(list_by_seg_id)

                # bother to coreg only if there is more than 1 segment
                if len(list_by_seg_id) > 1:

                    tmp_dir_coreg = create_tmp_dir_for_outfile(os.path.join(tmp_dir, 'coreg'))
                    # master dem is the first segment
                    dem_master = list_by_seg_id[0]
                    k = 0

                    # for csv summary: first segment
                    dem_master_id_split = os.path.basename(dem_master).split('_')
                    dem_master_id = '_'.join(dem_master_id_split[0:5])
                    dem_master_seg_nb = int(dem_master_id_split[5][3:])
                    writer.writerow(
                        [dem_master_id, dem_master_seg_nb, dem_master_seg_nb, 'first_seg', 'first_seg', 'first_seg',
                         'first_seg', str(k + 1)])

                    for i in range(len(list_by_seg_id) - 1):

                        dem_in = list_by_seg_id[i + 1]
                        print('Master segment is ' + dem_master)
                        print('Slave segment is ' + dem_in)

                        # for csv summary: all segments
                        dem_in_id_split = os.path.basename(dem_in).split('_')
                        dem_in_id = '_'.join(dem_in_id_split[0:5])
                        dem_in_seg_nb = int(dem_in_id_split[5][3:])
                        dem_master_id_split = os.path.basename(dem_master).split('_')
                        dem_master_seg_nb = int(dem_master_id_split[5][3:])
                        k0 = k

                        # deriving intersection
                        extent = inters_raster(dem_in, dem_master)

                        nb_init_pixel = valid_points_coreg(dem_master, dem_in)

                        # try coreg only if intersection is non void
                        if extent is not None and nb_init_pixel > nb_init_pixel_min:

                            tmp_dir_seg = create_tmp_dir_for_outfile(dem_in)

                            # do coreg, save results
                            try:
                                out_offs, final_rmse = dem_coregistration_custom(dem_master, dem_in,
                                                                                 outdir=tmp_dir_seg)[2:4]
                            except Exception:
                                out_offs = None
                                final_rmse = 999

                            if final_rmse < max_coreg_rmse:
                                str_out = 'good'
                            else:
                                str_out = 'breaks'

                            shutil.move(os.path.join(tmp_dir_seg, 'CoRegistration_Results.pdf'),
                                        os.path.join(dem_dir_out, 'coreg_results', str_out,
                                                     os.path.splitext(os.path.basename(dem_in))[0] + '_' +
                                                     os.path.splitext(os.path.basename(dem_master))[
                                                         0] + '_coreg_results.pdf'))
                            shutil.move(os.path.join(tmp_dir_seg, 'coreg_params.txt'),
                                        os.path.join(dem_dir_out, 'coreg_results', str_out,
                                                     os.path.splitext(os.path.basename(dem_in))[0] + '_' +
                                                     os.path.splitext(os.path.basename(dem_master))[
                                                         0] + '_coreg_params.txt'))

                            print('Final coregistration RMSE is:' + str(final_rmse))

                            # if rmse is good, keep doing mosaics with next segment: mosaic is new master
                            if final_rmse < max_coreg_rmse:
                                print('All good. Doing mosaic and iterating to next segment...')
                                tmp_mosaic_list = [dem_master, os.path.join(tmp_dir_seg,
                                                                            os.path.splitext(os.path.basename(dem_in))[
                                                                                0] + '_adj.tif')]

                                tmp_dem_coreg = os.path.join(tmp_dir_coreg, os.path.basename(list_by_seg_id[i + 1]))

                                dem_mosaic(tmp_mosaic_list, tmp_dem_coreg, erode_length=0, hole_fill_length=0,
                                           blending_method='mean')
                                dem_master = tmp_dem_coreg

                                # if this is the last segment, save mosaic and add to list of finalized mosaics
                                if i == (len(list_by_seg_id) - 2):
                                    k = k + 1
                                    final_strip = os.path.join(tmp_dir,
                                                               list_seg_id[list_list_by_seg_id.index(
                                                                   list_by_seg_id)] + '_mosaic' + str(
                                                                   k) + '_' + seg_suf)
                                    shutil.copyfile(dem_master, final_strip)
                                    list_final_strip.append(final_strip)

                            # otherwise, break the iterative mosaic-ing
                            else:
                                # if this is the last segment
                                if i == (len(list_by_seg_id) - 2):
                                    # if master is the previous segment alone
                                    if dem_master == list_by_seg_id[i]:
                                        print('RMSE is larger than limit of ' + str(
                                            max_coreg_rmse) + 'm: saving last segments independently...')
                                        list_final_strip.append(dem_in)
                                        list_final_strip.append(dem_master)
                                    # if master is a previous ongoing mosaic
                                    else:
                                        k = k + 1
                                        final_strip = os.path.join(tmp_dir,
                                                                   list_seg_id[list_list_by_seg_id.index(
                                                                       list_by_seg_id)] + '_mosaic' + str(
                                                                       k) + '_' + seg_suf)
                                        shutil.copyfile(dem_master, final_strip)
                                        print('RMSE is larger than limit of ' + str(
                                            max_coreg_rmse) + 'm: saving previous mosaic and last segment...')
                                        list_final_strip.append(dem_in)
                                        list_final_strip.append(final_strip)
                                # if there is still segments to go
                                else:
                                    # save current mosaic
                                    k = k + 1
                                    final_strip = os.path.join(tmp_dir,
                                                               list_seg_id[list_list_by_seg_id.index(
                                                                   list_by_seg_id)] + '_mosaic' + str(
                                                                   k) + '_' + seg_suf)
                                    shutil.copyfile(dem_master, final_strip)
                                    print('RMSE is larger than limit of ' + str(
                                        max_coreg_rmse) + 'm: saving previous mosaic, slave is new master...')
                                    list_final_strip.append(final_strip)
                                    # new master is the segment that did not mosaic
                                    dem_master = dem_in

                            remove_tmp_dir_for_outfile(dem_in)

                        # if no intersection and last segment
                        elif i == (len(list_by_seg_id) - 2):
                            out_offs = None
                            final_rmse = None
                            if dem_master == list_by_seg_id[i]:
                                print('No intersection between master and slave: saving last segments independently...')
                                list_final_strip.append(dem_in)
                                list_final_strip.append(dem_master)
                            # if master is a previous ongoing mosaic
                            else:
                                k = k + 1
                                final_strip = os.path.join(tmp_dir,
                                                           list_seg_id[list_list_by_seg_id.index(
                                                               list_by_seg_id)] + '_mosaic' + str(
                                                               k) + '_' + seg_suf)
                                shutil.copyfile(dem_master, final_strip)
                                print(
                                    'No intersection between master and slave: saving previous mosaic and last segment...')
                                list_final_strip.append(dem_in)
                                list_final_strip.append(final_strip)
                        # if no intersection and still segments to go
                        else:
                            out_offs = None
                            final_rmse = None
                            # save current mosaic
                            k = k + 1
                            final_strip = os.path.join(tmp_dir,
                                                       list_seg_id[list_list_by_seg_id.index(
                                                           list_by_seg_id)] + '_mosaic' + str(
                                                           k) + '_' + seg_suf)
                            shutil.copyfile(dem_master, final_strip)
                            print(
                                'No intersection between master and slave: saving previous mosaic, slave is new master...')
                            list_final_strip.append(final_strip)
                            # new master is the segment that did not mosaic
                            dem_master = dem_in

                        # write data in summary
                        if out_offs is None:
                            dx = 'no_inters'
                            dy = 'no_inters'
                            dz = 'no_inters'
                            rmse = 'no_inters'
                        else:
                            dx = str(out_offs[0])
                            dy = str(out_offs[1])
                            dz = str(out_offs[2])
                            rmse = str(final_rmse)

                        writer.writerow([dem_in_id, dem_in_seg_nb, dem_master_seg_nb, dx, dy, dz, rmse, str(k0 + 1)])

                    remove_tmp_dir_for_outfile(os.path.join(tmp_dir, 'coreg'))

                else:
                    list_final_strip.append(list_by_seg_id[0])

                    dem_master_id_split = os.path.basename(list_by_seg_id[0]).split('_')
                    dem_master_id = '_'.join(dem_master_id_split[0:5])
                    dem_master_seg_nb = int(dem_master_id_split[5][3:])
                    writer.writerow(
                        [dem_master_id, dem_master_seg_nb, dem_master_seg_nb, 'only_seg', 'only_seg', 'only_seg',
                         'only_seg', 'only_seg'])


    else:
        list_final_strip = list_tmp_dem

    # 3/ PROCESSING OF FINAL STRIPS

    final_dem = []
    for final_strip in list_final_strip:

        # 3.1/ FILTER STRIPS
        if filter_params is not None:
            sys.exit('No filter pre-defined for this DEM product.')

        # 3.2/ REPROJECT STRIPS

        # raw data is GeoTiff, 3413, 1 arc-sec and -9999 nodata_out
        if format_out == 'GTiff' and tgt_EPSG == 3031 and tgt_res is None and nodata_out is -9999:
            tmp_dem_proj = final_strip
        else:
            tmp_dem_proj = os.path.join(tmp_dir, os.path.splitext(os.path.basename(final_strip))[0] + '_proj.tif')
            warp_defaultUTM(final_strip, tmp_dem_proj, format_out, 3031, tgt_EPSG, tgt_res, nodata_out, interp_method)

        # 3.3/ ELLIPSOID OR GEOID STRIPS

        # raw data is ellipsoid WGS84
        if geoid:
            tmp_dem_geoid = os.path.join(tmp_dir, os.path.splitext(os.path.basename(final_strip))[0] + '_geoid.tif')
            ellipsoid_to_geoid(tmp_dem_proj, tmp_dem_geoid)
        else:
            tmp_dem_geoid = tmp_dem_proj

        final_dem.append(tmp_dem_geoid)

    for dem in final_dem:
        shutil.copyfile(dem, os.path.join(dem_dir_out, os.path.basename(dem)))

    remove_tmp_dir_for_outfile(os.path.join(dem_dir_out, 'all_strips'))

    print('Fin.')