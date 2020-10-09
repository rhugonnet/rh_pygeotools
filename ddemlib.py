# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 11:55:37 2018

@author: hugonnet

DDEM LIBRARY:
Library of Python functions for manipulating DEM differences
"""
import os, sys, shutil
import csv
import numpy as np
import pandas as pd
from numpy.polynomial.polynomial import polyfit, polyval
import random
from vectlib import simplify_shp_fn, buffer_shp_fn, inters_shp_fn, union_shp_fn, poi_polygon, isempty_firstfeat, copy_shp_fn, extent_shp_ref
from rastlib import rasterize_shp, proximity_shp, polygonize_fn, write_nanarray, read_nanarray, proximity_rast_fn, pixel_size, create_mem_raster_on_ref
from shlib import create_tmp_dir_for_outfile, remove_tmp_dir_for_outfile
from scipy.ndimage.filters import convolve
import scipy.stats as st
from fillalglib import floodfill_discontinuous
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
from subprocess import Popen
from pybob.GeoImg import GeoImg
from pybob.plot_tools import plot_ddem_results, plot_polygon_df
import geopandas as gpd

def ddem_discrete_hypso(ddem,dem,mask,gsd,proxi=None,bin_type='fixed',bin_val=50.,filt='5NMAD'):

    final_mask = np.logical_and(np.logical_and(np.isfinite(ddem), np.isfinite(dem)),mask)

    dem_on_mask = dem[final_mask]
    ddem_on_mask = ddem[final_mask]
    if proxi is not None:
        proxi_on_mask = proxi[final_mask]

    ddem_out = np.copy(ddem)

    min_elev = np.min(dem_on_mask) - (np.min(dem_on_mask) % bin_val)
    max_elev = np.max(dem_on_mask) + 1

    if bin_type == 'fixed':
        bin_final = bin_val
    elif bin_type == 'percentage':
        bin_final = np.ceil(bin_val / 100. * (max_elev - min_elev))
    else:
        sys.exit('Bin type not recognized.')

    bins_on_mask = np.arange(min_elev, max_elev, bin_final)

    nb_bin = len(bins_on_mask)

    elev_bin = np.zeros(nb_bin)*np.nan
    nmad_bin = np.zeros(nb_bin)*np.nan
    med_bin = np.zeros(nb_bin)*np.nan
    std_bin = np.zeros(nb_bin)*np.nan
    area_tot_bin = np.zeros(nb_bin)*np.nan
    area_meas_bin = np.zeros(nb_bin)*np.nan
    prox = np.zeros(nb_bin)*np.nan

    for i in np.arange(nb_bin):

        idx_bin = np.array(dem_on_mask >= bins_on_mask[i]) & np.array(
            dem_on_mask < (bins_on_mask[i] + bin_final))
        idx_orig = np.array(dem >= bins_on_mask[i]) & np.array(
            dem < (bins_on_mask[i] + bin_final)) & mask
        area_tot_bin[i] = np.count_nonzero(idx_orig)*gsd**2
        area_meas_bin[i] = np.count_nonzero(idx_bin)*gsd**2
        elev_bin[i] = bins_on_mask[i] + bin_final / 2.
        dh_bin = ddem_on_mask[idx_bin]
        if proxi is not None:
            proxi_bin = proxi_on_mask[idx_bin]
            if len(proxi_bin[~np.isnan(proxi_bin)])>0:
                prox[i] = np.nanmax(proxi_bin)

        if len(dh_bin[~np.isnan(dh_bin)]) > 0:

            std_bin[i] = np.nanstd(dh_bin)
            med_bin[i] = np.nanmedian(dh_bin)
            if filt:
                median_temp = np.nanmedian(dh_bin)
                MAD_temp = np.nanmedian(np.absolute(dh_bin[~np.isnan(dh_bin)] - median_temp))
                NMAD_temp = 1.4826 * MAD_temp
                nmad_bin[i] = NMAD_temp
                # dh_bin[np.absolute(dh_bin - median_temp) > 5 * NMAD_temp] = np.NaN

                ddem_out[idx_orig & np.array(np.absolute(ddem_out - median_temp) > 5 * NMAD_temp)] = np.nan

    return ddem_out, elev_bin, med_bin, std_bin, nmad_bin, area_tot_bin, area_meas_bin, prox


def plot_hypso_fit(ddem_masked,dem_masked,elev,med,nmad,std,elevfit,poly_order,pp):

    mykeep = np.logical_and(np.isfinite(ddem_masked),np.isfinite(dem_masked))
    H = dem_masked[mykeep]
    dH = ddem_masked[mykeep]
    sampsize = 25000
    if H.size > sampsize:
        mysamp = np.random.randint(0, H.size, sampsize)
    else:
        mysamp = np.arange(0, H.size)

    newelev=np.arange(min(elev),max(elev),1)
    interp_dH = polyval(newelev,elevfit)
    first_der=np.polyder(elevfit)
    second_der=np.polyder(first_der)
    der1_dH = polyval(newelev,first_der)
    der2_dH = polyval(newelev,second_der)

    fig = plt.figure(figsize=(7, 5), dpi=300)
    plt.title('Hypsometric distribution of elevation change', fontsize=14)
    plt.plot(H[mysamp], dH[mysamp], '^', ms=0.75, color='0.5', rasterized=True, fillstyle='full',
             label="Raw [samples]")
    plt.plot(elev, med, '-', ms=2, color='0.15', label='median')

    plt.plot(elev, nmad, 'r--', ms=2, label="nmad")
    plt.plot(elev, std, 'm--', ms=2, label="std")
    plt.plot(newelev,interp_dH,'b--',ms=2,label='polynomial fit: order '+str(poly_order))

    plt.xlim(np.min(elev), np.max(elev))
    plt.ylim(np.min(med) - np.max(std),np.max(std))

    plt.xlabel(' Elevation [m]')
    plt.ylabel('dH [m]')
    plt.legend(loc='lower right')

    pp.savefig(fig, dpi=300)

    fig2 = plt.figure(figsize=(7, 5), dpi=300)
    plt.title('Hypsometric distribution of elevation change', fontsize=14)
    plt.plot(newelev,der1_dH,'g--',ms=2,label='first derivative')
    plt.xlim(np.min(elev), np.max(elev))
    plt.ylim(np.min(der1_dH),np.max(der1_dH))
    plt.xlabel(' Elevation [m]')
    plt.ylabel('dH [m]')
    plt.legend(loc='lower right')
    pp.savefig(fig2, dpi=300)

    fig3 = plt.figure(figsize=(7, 5), dpi=300)
    plt.title('Hypsometric distribution of elevation change', fontsize=14)
    plt.plot(newelev,der2_dH,'c--',ms=2,label='second derivative')
    plt.xlim(np.min(elev), np.max(elev))
    plt.ylim(np.min(der2_dH),np.max(der2_dH))
    plt.xlabel(' Elevation [m]')
    plt.ylabel('dH [m]')
    plt.legend(loc='lower right')
    pp.savefig(fig3, dpi=300)


def ddem_poly_hypso(ddem,dem,mask,gsd,poly_order=5,pp=None,filter_3nmad=True,get_gap_filled=False,get_elev_residual=False):

    ddem_filtered, elev, med, std, nmad, area_tot, area_meas = ddem_discrete_hypso(ddem, dem, mask,gsd=gsd,bin_val=10.)

    if filter_3nmad:
        final_mask = np.logical_and(np.logical_and(np.isfinite(ddem_filtered), np.isfinite(dem)),mask)
    else:
        final_mask = np.logical_and(np.logical_and(np.isfinite(ddem), np.isfinite(dem)), mask)

    perc_mask=10.
    critic_elev=np.nanmin(dem[mask])+perc_mask*(np.nanmax(dem[mask])-np.nanmin(dem[mask]))/100

    final_mask = np.logical_and(final_mask,(dem>critic_elev))

    ddem_data = ddem[final_mask]
    dem_data = dem[final_mask]

    elevfit = polyfit(elev[elev>critic_elev], med[elev>critic_elev], poly_order)
    hole_inds = np.where(np.logical_and(np.logical_and(np.isnan(ddem),
                                 np.isfinite(dem)),mask))

    ddem_out = np.copy(ddem)

    if get_gap_filled:
        hole_elevs = dem[hole_inds]
        interp_points = polyval(hole_elevs, elevfit)
        ddem_out[hole_inds] = interp_points

    if get_elev_residual:
        #get polynomial fit on nmad elevation trend
        nmadfit = polyfit(elev,std,3)

        apply_mask = np.logical_and(np.logical_and(np.isfinite(ddem), np.isfinite(dem)), mask)
        dem_apply = dem[apply_mask]
        # difference each pixel elevation with polynomial fit
        interp_elev = polyval(dem_apply, elevfit)
        ddem_out[apply_mask] -= interp_elev
        # standardized by elevation variance (nmad)
        interp_nmad = polyval(dem_apply, nmadfit)
        ddem_out[apply_mask] /= interp_nmad

    if pp is not None:
        plot_hypso_fit(ddem_data,dem_data,elev,med,nmad,std,elevfit,poly_order,pp)

    return ddem_out, elev, med, std, nmad, elevfit


def plot_med_hypso_fit(ddem_masked,dem_masked,elev,med,nmad,std,pp):

    mykeep = np.logical_and(np.isfinite(ddem_masked),np.isfinite(dem_masked))
    H = dem_masked[mykeep]
    dH = ddem_masked[mykeep]
    sampsize = 25000
    if H.size > sampsize:
        mysamp = np.random.randint(0, H.size, sampsize)
    else:
        mysamp = np.arange(0, H.size)

    #temporary save sampled points for MdG figure
    # df=pd.DataFrame()
    # df=df.assign(rdn_elev=H[mysamp],rdn_dh=dH[mysamp])
    # df.to_csv('/home/atom/ongoing/std_err/data_vhr/etienne_mb/geophys_var/df_rand_dh.csv')

    fig = plt.figure(figsize=(7, 5), dpi=300)
    plt.title('Hypsometric distribution of elevation change', fontsize=14)
    plt.plot(H[mysamp], dH[mysamp], '^', ms=0.75, color='0.5', rasterized=True, fillstyle='full',
             label="Raw [samples]")
    plt.plot(elev, med, '-', ms=2, color='0.15', label='median')

    plt.plot(elev, nmad, 'r--', ms=2, label="nmad")
    plt.plot(elev, std, 'm--', ms=2, label="std")

    plt.xlim(np.nanmin(elev), np.nanmax(elev))
    plt.ylim(np.nanmin(med) - np.nanmax(std),np.nanmax(std))

    plt.xlabel(' Elevation [m]')
    plt.ylabel('dH [m]')
    plt.legend(loc='lower right')

    pp.savefig(fig, dpi=300)

def ddem_med_hypso(ddem,dem,mask,gsd,pp=None,proxi=None,filter_3nmad=True,get_gap_filled=False,get_elev_residual=False):

    ddem_filtered, elev, med, std, nmad, area_tot, area_meas, prox = ddem_discrete_hypso(ddem, dem, mask,gsd=gsd,proxi=proxi,bin_val=50.)

    if filter_3nmad:
        final_mask = np.logical_and(np.logical_and(np.isfinite(ddem_filtered), np.isfinite(dem)),mask)
    else:
        final_mask = np.logical_and(np.logical_and(np.isfinite(ddem), np.isfinite(dem)), mask)

    # perc_mask=10.
    # critic_elev=np.nanmin(dem[mask])+perc_mask*(np.nanmax(dem[mask])-np.nanmin(dem[mask]))/100
    #
    # final_mask = np.logical_and(final_mask,(dem>critic_elev))

    ddem_data = ddem[final_mask]
    dem_data = dem[final_mask]

    res = np.copy(ddem)
    res_stdized = np.copy(ddem)

    # if get_gap_filled:
    #     hole_elevs = dem[hole_inds]
    #     interp_points = polyval(hole_elevs, elevfit)
    #     ddem_out[hole_inds] = interp_points

    if get_elev_residual:

        apply_mask = np.logical_and(np.logical_and(np.isfinite(ddem), np.isfinite(dem)), mask)
        dem_apply = dem[apply_mask]
        # difference each pixel elevation with median interp
        interp_elev = np.interp(dem_apply, elev, med)
        res[apply_mask] -= interp_elev
        # standardized by elevation variance (nmad)
        interp_nmad = np.interp(dem_apply, elev, std)
        res_stdized = res
        res_stdized[apply_mask] /= interp_nmad

    if pp is not None:
        plot_med_hypso_fit(ddem_data,dem_data,elev,med,nmad,std,pp)

    return res, res_stdized, elev, med, std, nmad, area_tot, area_meas, prox


def get_geophys_var_hypso(fn_ddem,fn_dem,fn_shp,out_dir,path_to_r_geophys):


    pp = PdfPages(os.path.join(out_dir, 'hypsometric_fit_results.pdf'))

    ddem=read_nanarray(fn_ddem)
    ddem[np.absolute(ddem)>60]=np.nan
    # ddem = ddem*12.
    dem=read_nanarray(fn_dem)
    mask = rasterize_shp(fn_shp, fn_dem)
    gsd = pixel_size(fn_ddem)
    fn_proxi=os.path.join(out_dir,'proxi.tif')
    proxi = proximity_shp(fn_shp,fn_ddem,type_prox='interior')

    #first get residuals of poly fit
    res, res_stdized, elev, med, std, nmad, area_tot, area_meas, prox = ddem_med_hypso(ddem,dem,mask,gsd,pp=pp,proxi=proxi,get_elev_residual=True)
    plt.close('all')

    fn_mask=os.path.join(out_dir,'mask.tif')
    write_nanarray(fn_mask,fn_ddem,mask)
    fn_res=os.path.join(out_dir,'residual.tif')
    fn_res_stdized = os.path.join(out_dir,'residual_standardized.tif')
    write_nanarray(fn_res,fn_ddem,res)
    write_nanarray(fn_res_stdized,fn_ddem,res_stdized)

    mask_geo = GeoImg(fn_mask)
    res_geo = GeoImg(fn_res)
    res_stdized_geo = GeoImg(fn_res_stdized)
    ddem_geo = GeoImg(fn_ddem)
    dem_geo=GeoImg(fn_dem)
    # res_geo.img[np.invert(mask)] = np.nan
    extent = extent_shp_ref(fn_shp,fn_dem)
    crop_res = res_geo.crop_to_extent([extent[0],extent[2],extent[1],extent[3]])
    crop_res_stdized = res_stdized_geo.crop_to_extent([extent[0],extent[2],extent[1],extent[3]])
    crop_ddem = ddem_geo.crop_to_extent([extent[0],extent[2],extent[1],extent[3]])
    crop_mask = mask_geo.crop_to_extent([extent[0],extent[2],extent[1],extent[3]])
    crop_dem = dem_geo.crop_to_extent([extent[0],extent[2],extent[1],extent[3]])

    fn_crop_res_stdized = os.path.join(out_dir,'res_stdized_cropped.tif')
    fn_crop_mask = os.path.join(out_dir,'mask_cropped.tif')
    fn_crop_dem = os.path.join(out_dir,'dem_cropped.tif')

    crop_res_stdized.img[crop_mask.img != 1]=np.nan
    crop_res_stdized.write(fn_crop_res_stdized)
    crop_mask.write(fn_crop_mask)
    crop_dem.write(fn_crop_dem)

    crop_res.img[crop_mask.img != 1]=np.nan
    crop_res_stdized.img[crop_mask.img != 1]=np.nan
    # crop_ddem.img = 12*crop_ddem.img


    clim_ddem_raw = np.nanmax(np.absolute(med))

    outline_gla = gpd.read_file(fn_shp)
    fig, _ = plot_polygon_df(outline_gla, edgecolor='k', lw=2, alpha=0.5)
    plt.title('Outline')
    pp.savefig(fig, dpi=300)

    fig = plot_ddem_results(crop_ddem, clim=(-clim_ddem_raw,clim_ddem_raw), colormap='Spectral')[0]
    plt.title('Elevation change [m] (Large scale)')
    pp.savefig(fig, dpi=300)

    fig = plot_ddem_results(crop_ddem, clim=(-3, 3), colormap='Spectral')[0]
    plt.title('Elevation change [m] (Thin scale)')
    pp.savefig(fig, dpi=300)

    clim_res = np.nanmean(np.absolute(nmad))

    fig = plot_ddem_results(crop_res, clim=(-clim_res,clim_res),colormap='Spectral')[0]
    plt.title('Hypsometric residual of elevation change [m] \n (Elevation change minus hypsometric median)')
    pp.savefig(fig, dpi=300)

    fig = plot_ddem_results(crop_res_stdized, clim=(-1, 1), colormap='Spectral')[0]
    plt.title('Standardized hypsometric residual of elevation change [no unit] \n (Elevation change minus hypsometric median divided by hypsometric nmad)')
    pp.savefig(fig, dpi=300)

    pp.close()
    plt.close('all')
    os.remove(fn_res)
    os.remove(fn_mask)
    os.remove(fn_res_stdized)


    #normalize elevation
    max_elev=np.nanmax(elev)
    min_elev=np.nanmin(elev)
    elev_n = (elev - min_elev) / (max_elev - min_elev)

    #normalize dh
    max_dh=np.nanmax(med)
    min_dh=np.nanmin(med)
    accu_elev = min_elev + 80*(max_elev - min_elev)/100
    tmp_max_dh=np.nanmean(med[elev>accu_elev]) #use mean of accumulation instead of max
    if np.abs((np.nanmax(med)-tmp_max_dh)/(max_dh-min_dh))<0.3:
        max_dh = tmp_max_dh
    med_n = (min_dh - med) / (max_dh - min_dh)
    std_n = std / (max_dh - min_dh)
    nmad_n = nmad / (max_dh - min_dh)

    #write normalized data
    elev_rs = np.arange(0,1,0.01)
    med_rs = np.interp(elev_rs,elev_n,med_n)
    std_rs = np.interp(elev_rs,elev_n,std_n)
    nmad_rs = np.interp(elev_rs,elev_n,nmad_n)
    area_rs = np.interp(elev_rs,elev_n,area_tot)
    df = pd.DataFrame()
    df= df.assign(norm_elev=elev_rs,norm_med_dh=med_rs,norm_std_dh=std_rs,norm_nmad_rs=nmad_rs,area_rs=area_rs)
    df_raw = pd.DataFrame()
    df_raw =df_raw.assign(elev=elev,med_dh=med,std_dh=std,nmad_dh=nmad,area_tot=area_tot,area_meas=area_meas,prox=prox)

    df.to_csv(os.path.join(out_dir,'df_norm_dh_elev.csv'))
    df_raw.to_csv(os.path.join(out_dir,'df_raw_dh_elev.csv'))

    ddem = dem = mask = res = res_stdized = crop_mask = crop_res_stdized = crop_res = crop_ddem = crop_dem = ddem_geo = dem_geo = res_geo = res_stdized_geo = None

    #get variogram with moving elevation window from R
    # cmd = 'Rscript '+path_to_r_geophys+' -d '+fn_crop_dem+' -r '+fn_crop_res_stdized+' -m '+fn_crop_mask+' -v Exp -o '+out_dir
    # fn_log = os.path.join(out_dir,'r_geophys.log')
    # log=open(fn_log,'w')
    # p=Popen(cmd,stdout=log,stderr=log,shell=True)
    # p.wait()
    # log.close()

    os.remove(fn_crop_dem)
    os.remove(fn_crop_res_stdized)
    os.remove(fn_crop_mask)


def kernel_tvol_outline_adjust(fn_ddem,fn_dem,fn_outline_in,fn_outline_out,fn_mask_stable,fn_mask_gla,tolerance_factor=5.,kernel_size=9,only_thinning=True):

    def gkern(kernlen=21):

        #source: https://stackoverflow.com/questions/29731726/how-to-calculate-a-gaussian-kernel-matrix-efficiently-in-numpy
        """Returns a 2D Gaussian kernel."""

        lim = kernlen // 2 + (kernlen % 2) / 2
        x = np.linspace(-lim, lim, kernlen + 1)
        kern1d = np.diff(st.norm.cdf(x))
        kern2d = np.outer(kern1d, kern1d)
        return kern2d / kern2d.sum()

    #inside parameters based on prior knowledge:
    tmp_dir=create_tmp_dir_for_outfile(fn_outline_out)

    #read
    ddem = read_nanarray(fn_ddem)
    ddem[np.absolute(ddem)>60]=np.nan
    dem = read_nanarray(fn_dem)
    mask_other_gla = (read_nanarray(fn_mask_stable) == 1)
    mask_gla = (read_nanarray(fn_mask_gla) == 1)
    mask_stable = np.invert(np.logical_or(mask_gla,mask_other_gla))
    res = pixel_size(fn_dem)

    elev, med, std, nmad = ddem_discrete_hypso(ddem, dem, mask_gla,gsd=res,bin_val=10)[1:5]

    #first, keep stable terrain
    keep_stable = np.logical_and(mask_stable,np.isfinite(ddem))
    ddem_on_stable=ddem[keep_stable]

    #get nmad as representative statistic of stable terrain
    median_stable = np.nanmedian(ddem_on_stable)
    mad_stable = np.nanmedian(np.absolute(ddem_on_stable[~np.isnan(ddem_on_stable)] - median_stable))
    nmad_stable = 1.4826 * mad_stable

    # #filter heavy outliers to get a fair representation of overall stable terrain (without unmapped glaciers, landslides, etc...)
    # ddem_on_stable[np.array(np.absolute(ddem_on_stable - median_stable) > 3 * nmad_stable)] = np.nan

    #rasterize outline
    # mask_outline = rasterize_shp(fn_outline_in,fn_ddem)

    #get an idea of thinning rate and amplitude of adjustment
    # dh_on_outline = ddem_poly_elev(ddem,dem,5,mask_outline)

    outline_timelapse = 10
    ddem_timelapse=10
    # detect elev above which thinning is not significant anymore so we don't care
    critic_elev = np.nanmax(elev[np.absolute(med)+3*np.absolute(nmad)>tolerance_factor*nmad_stable])
    filt_elev = (elev<=critic_elev)
    # remove decreasing dh signal around the end of ablation area to assess effective dh with elevation
    filt_diff_dh = (np.diff(med) > 0)
    filt_diff_dh=np.append(filt_diff_dh,True)
    filt=np.logical_and(filt_elev,filt_diff_dh)
    new_elev=elev[filt]
    new_med=med[filt]
    #range of proximity value to search in as a function of elevation
    critic_proximity = np.absolute(new_med)*outline_timelapse*ddem_timelapse #TODO: ideally, multiply that by a timelapse between outline and first DEM, make sure ddem is in m and not m.y-1, and then apply a factor with slope convolution in the proximity area
    interp_critic_proximity = np.reshape(np.interp(dem.flatten(),new_elev,critic_proximity,right=0,left=np.nanmax(critic_proximity)),np.shape(dem)) #TODO: could extrapolate linearly with elevation to assess lower dh?

    #define mask of areas for possible outline adjustment
    # mask_surround_elev = (dem>min_elev-mvp_elev_value) & (dem<=critic_elev)
    proximity_rast = proximity_shp(fn_outline_in,fn_ddem)
    mask_proximity = proximity_rast<interp_critic_proximity

    #define mask for neighboring glaciers
    tmp_rast_othergla = os.path.join(tmp_dir, 'tmp_rast_othergla.tif')
    tmp_prox_othergla = os.path.join(tmp_dir, 'tmp_prox_othergla.tif')
    rast_othergla = np.ones(np.shape(mask_other_gla))
    rast_othergla[~mask_other_gla] = 0
    write_nanarray(tmp_rast_othergla, fn_dem, rast_othergla)
    proximity_rast_fn(tmp_rast_othergla, tmp_prox_othergla)
    prox_othergla = read_nanarray(tmp_prox_othergla)

    max_prox_othergla = 4 * kernel_size * res
    mask_max_prox_othergla = prox_othergla > max_prox_othergla

    mask_maybe = np.logical_and(mask_proximity, mask_max_prox_othergla)

    #get upper elevation/center of glacier by using intersect with original outline with buffer/polygonize excluding proximity mask
    tmp_buffered = os.path.join(tmp_dir,'tmp_buff.shp')
    buffer_shp_fn(fn_outline_in,tmp_buffered,3*res)
    buff_rast=rasterize_shp(tmp_buffered,fn_ddem)
    buff_rast[mask_maybe]=0
    buff_mask=np.ones(np.shape(buff_rast))
    buff_mask[~buff_rast]=-9999
    tmp_buffered_rast= os.path.join(tmp_dir,'tmp_buff_mask.tif')
    tmp_upper_outline=os.path.join(tmp_dir,'tmp_upper.shp')
    write_nanarray(tmp_buffered_rast,fn_ddem,buff_mask)
    #if there is an area untouched that is kept to union with adjusted polygon at the end
    nb_pixel_kept=len(buff_mask[buff_mask==1])
    if nb_pixel_kept>0:
        polygonize_fn(tmp_buffered_rast,tmp_upper_outline)
        tmp_inters=os.path.join(tmp_dir,'tmp_inters.shp')
        # inters_shp_fn(fn_outline_in,tmp_upper_outline,tmp_inters)

        #get a common area to merge with final polygon by buffering previous outline
        tmp_upper_buff=os.path.join(tmp_dir,'tmp_buff_upper.shp')
        tmp_upper_buff_inters=os.path.join(tmp_dir,'tmp_buff_upper_inters.shp')
        buffer_shp_fn(tmp_upper_outline,tmp_upper_buff,2*kernel_size*res)
        inters_shp_fn(fn_outline_in,tmp_upper_buff,tmp_upper_buff_inters)
        mask_upper_buff = rasterize_shp(tmp_upper_buff_inters,fn_dem)
    #if there is no area untouched, find pole of inaccessibility
    else:
        tmp_poi=os.path.join(tmp_dir,'poi.shp')
        poi_polygon(fn_outline_in,tmp_poi,res)
        mask_upper_buff = rasterize_shp(tmp_poi,fn_dem,all_touched=True)

    ddem_maybe=ddem
    ddem_maybe[~mask_maybe] = np.nan
    if only_thinning:
        ddem_maybe[ddem_maybe>5*nmad_stable] = np.nan

    #use kernel to ensure volume change at a given pixel is significant
    conv = convolve(ddem_maybe,gkern(kernel_size),mode='nearest')
    if only_thinning:
        mask_tvol = (conv < -tolerance_factor*nmad_stable)
    else:
        mask_tvol = (np.absolute(conv) > tolerance_factor*nmad_stable)

    #use filling algorithm from common mask to get only spatially continuous adjusted area with interior of glacier
    mask_common = np.logical_and(mask_upper_buff, mask_tvol)
    mask_flood = floodfill_discontinuous(mask_tvol,mask_start=mask_common)
    mask_flood = mask_flood.astype('float')
    mask_flood[mask_flood==0]=-9999

    fn_mask_flood =os.path.join(tmp_dir,'mask_flood.tif')
    fn_shp_flood = os.path.join(tmp_dir,'shp_flood.shp')
    fn_smooth_flood = os.path.join(tmp_dir,'shp_flood_smooth.shp')

    write_nanarray(fn_mask_flood,fn_ddem,mask_flood)

    #polygonize final mask
    polygonize_fn(fn_mask_flood,fn_shp_flood)

    #smooth polygon
    simplify_shp_fn(fn_shp_flood,fn_smooth_flood,tolerance=res)

    #union with inside
    if nb_pixel_kept>0 and not isempty_firstfeat(fn_smooth_flood):
        union_shp_fn(fn_smooth_flood,tmp_upper_buff_inters,fn_outline_out)
    elif not isempty_firstfeat(fn_smooth_flood):
        copy_shp_fn(fn_smooth_flood,fn_outline_out)
    else:
        #TODO: actually need to change projection in this step or it stays that of the input shapefile
        copy_shp_fn(fn_outline_in,fn_outline_out)

    # remove_tmp_dir_for_outfile(fn_outline_out)


def simulate_sin_bias(fn_dem_ref,fn_along_bias,along_angle,ampli,freq):

    ds_along =create_mem_raster_on_ref(fn_dem_ref)
    res = pixel_size(fn_dem_ref)
    bin_size = 3*res

    dh_along = ds_along.GetRasterBand(1).ReadAsArray()
    nx, ny = np.shape(dh_along)
    xx, yy = np.meshgrid(np.arange(ny), np.arange(nx))
    altrack = -xx * np.sin(np.deg2rad(along_angle)) + yy * np.cos(np.deg2rad(along_angle))

    #correcting bias

    # filt = np.isfinite(dh_along)
    # xd = altrack[filt].flatten()
    # yd = dh_along[filt].flatten()
    # df = pd.DataFrame({'dh': yd, 'dist_altrack': xd})
    # df['dist_altrack'] = pd.cut(df['dist_altrack'], bins=np.arange(min(xd), max(xd) + bin_size, bin_size),
    #                             labels=np.arange(min(xd) + 50. / 2, max(xd) + bin_size / 2, bin_size))
    # df2 = df.groupby(['dist_altrack']).mean()
    # xtest = df2.index.values
    # ytest = df2['dh']

    #can occur when masking
    # filt_missing = np.logical_and(np.isfinite(xtest),np.isfinite(ytest))
    # keepxtest = xtest[filt_missing]
    # keepytest = ytest[filt_missing]

    # coeffs_al = np.polyfit(keepxtest, keepytest, 5)
    # poly_altrack = np.poly1d(coeffs_al)  # define a polynomial function
    # dh_corr = dh_along - poly_altrack(altrack)

    #simulating bias

    dh_simu =  ampli*np.sin(altrack*res*2*np.pi/freq)

    write_nanarray(fn_along_bias,fn_dem_ref,dh_simu)

def patches_method(fn_ddem,fn_mask_stable,perc_min_valid,averaging_area,method='circular',nmax=1000.):

    def create_circular_mask(h, w, center=None, radius=None):

        if center is None:  # use the middle of the image
            center = [int(w / 2), int(h / 2)]
        if radius is None:  # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w - center[0], h - center[1])

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

        mask = dist_from_center <= radius
        return mask

    ddem = read_nanarray(fn_ddem)
    res = pixel_size(fn_ddem)
    mask_stable = (read_nanarray(fn_mask_stable) > 0)

    valid_mask = np.logical_and(np.isfinite(ddem),mask_stable)

    ddem[~valid_mask] = np.nan
    nx, ny = np.shape(ddem)
    count = len(ddem[~np.isnan(ddem)])
    print('Number of valid pixels: ' +str(count))
    nb_cadrant = int(np.floor(np.sqrt((count * res * res) / averaging_area) + 1))
    #rectangular
    nx_sub = int(np.floor((nx - 1) / nb_cadrant))
    ny_sub = int(np.floor((ny - 1) / nb_cadrant))
    #radius size for a circular patch
    rad = int(np.floor(np.sqrt(averaging_area/res*res*np.pi)))

    tile = []
    mean_patch = []
    med_patch = []
    std_patch = []
    nb_patch = []

    list_cadrant = [[i,j] for i in range(nb_cadrant) for j in range(nb_cadrant)]
    u=0
    while len(list_cadrant)>0 and u<nmax:

        check = 0
        while check == 0:
            rand_cadrant = random.randint(0, len(list_cadrant)-1)

            i = list_cadrant[rand_cadrant][0]
            j = list_cadrant[rand_cadrant][1]

            check_x = int(np.floor(nx_sub*(i+1/2)))
            check_y = int(np.floor(ny_sub*(j+1/2)))
            if mask_stable[check_x,check_y]:
                check = 1

        list_cadrant.remove(list_cadrant[rand_cadrant])

        tile.append(int(str(i) + str(j)))
        if method == 'rectangular':
            patch = ddem[nx_sub * i:nx_sub * (i + 1), ny_sub * j:ny_sub * (j + 1)].flatten()
        elif method == 'circular':
            center_x = np.floor(nx_sub*(i+1/2))
            center_y = np.floor(ny_sub*(j+1/2))
            mask = create_circular_mask(nx,ny,center=[center_x,center_y],radius=rad)
            # patch = ddem[center_x-rad:center_x+rad,center_y-rad:center_y+rad]
            patch = ddem[mask]
        else:
            print('Patch method must be rectangular or circular.')
            sys.exit()

        nb_pixel_total = len(patch)
        nb_pixel_valid = len(patch[np.isfinite(patch)])
        if nb_pixel_valid > np.ceil(perc_min_valid / 100. * nb_pixel_total):
            print('Here : '+str(u))
            u=u+1
            mean_patch.append(np.nanmean(patch))
            med_patch.append(np.nanmedian(patch))
            std_patch.append(np.nanstd(patch))
            nb_patch.append(nb_pixel_valid)

    tile = np.array(tile)
    mean_patch = np.array(mean_patch)
    med_patch = np.array(med_patch)
    std_patch = np.array(std_patch)
    nb_patch = np.array(nb_patch)

    return tile, mean_patch, med_patch, std_patch, nb_patch

def patches_areas_to_csv(fn_ddem,fn_mask_stable,fn_csv_out,perc_min_valid,list_areas,method='circular',nmax=1000):

    with open(fn_csv_out, 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow([str(x) for x in list_areas])
        list_err = []
        list_nmad = []
        list_mean = []
        list_nb = []
        for averaging_area in list_areas:
            print('Working on area: ' +str(averaging_area))
            tile, mean, med, std, nb = patches_method(fn_ddem, fn_mask_stable,perc_min_valid=perc_min_valid,averaging_area=averaging_area,method=method,nmax=nmax)
            med_area = np.nanmedian(mean)
            mad_area = np.nanmedian(np.absolute(mean[~np.isnan(mean)] - med_area))
            nmad_area = 1.4826 * mad_area
            list_err.append(np.nanstd(mean))
            list_nmad.append(nmad_area)
            list_mean.append(np.nanmean(mean))
            list_nb.append(len(mean[~np.isnan(mean)]))

        print(list_err)
        writer.writerow([str(x) for x in list_err])
        writer.writerow([str(x) for x in list_nmad])
        writer.writerow([str(x) for x in list_mean])
        writer.writerow([str(x) for x in list_nb])


# if __name__ == '__main__':
#     fn_dem = '/home/atom/ongoing/RGI11/dh_RGI6.0_2000.5-2016.8/Trend_dh_RGI6.0_NotUsing_N46E008_weighted_2000.5-2016.8_BadTrendExcluded_merged_warped.tif'
#     fn_along_bias = '/home/atom/ongoing/test.tif'
#     along_angle = 23
#     ampli=1
#     freq=5000


