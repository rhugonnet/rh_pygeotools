

#in construction
from __future__ import print_function
import sys, os
import numpy as np
from math import pi, sqrt, exp
from statlib import lowess_homemade_kern, kernel_exclass, kernel_exp, kernel_gaussian
import pandas as pd

def old_mb_hypso_3NMAD(dh_trend, DEM, mask, ice_dens, min_alt, max_alt, pixel_size, bin_alti = 100, Id = 0):
    """
    Compute the mass balance of glacier defined by "mask" (0 = OFF GLA, any other value means GLA) and elevation changes "dh_trend".
    A mean elevation change is computed for each altitude band calculated from "DEM"
    Inputs :
    - dh_trend : array, map of elevation changes
    - mask : array, mask of the glacier, of same size as dh, with 0 outside of the glacier, >0 on the glacier
    - DEM : array, same size as dh, altitude of each point
    - pixel_size : f, pixel size in x/y direction in km, default = 0.03 (ASTER)
    - bin_alti : int, width of the elevation band, default = 100 m
    - Id : int, ID of the glacier for which you want to calculate the MB, default = 0 (glacier wide)

    Outputs:
    - np.nanmean(dh_trend_ON_GLA) = MB calculated after replacing all the outliers for GLA terrain
    - OUTPUT = table designed to match Etienne's template
    - dh_trend_ON_GLA = array, nan OFF GLA and filled with the mean per altitude band
    """

    mask_ON_GLA = np.ones(np.shape(mask))
    if Id == 0:
        mask_ON_GLA[mask != 1] = np.nan # 0 means OFF GLA
    else:
        mask_ON_GLA[mask != Id] = np.nan # keeps only

#    mask_OFF_GLA = np.ones(np.shape(mask))
#    mask_OFF_GLA[mask > 0.] = np.nan

    dh_trend_ON_GLA = dh_trend*mask_ON_GLA
    DEM_ON_GLA = DEM*mask_ON_GLA

    if min_alt is None:
        min_alt = np.nanmin(DEM_ON_GLA) - 0.5 * bin_alti
        max_alt = np.nanmax(DEM_ON_GLA) + 0.5 * bin_alti

    print(min_alt)
    print(max_alt)
    print(len(DEM_ON_GLA[~np.isnan(DEM_ON_GLA)]))

    bins_ON_GLA = np.arange(min_alt, max_alt, bin_alti)

    Mean_alti = np.zeros(len(bins_ON_GLA))*np.nan
    Mediane_Dh_all = np.zeros(len(bins_ON_GLA))*np.nan
    Mean_Dh_3nmad = np.zeros(len(bins_ON_GLA))*np.nan
    Std_dh_3nmad = np.zeros(len(bins_ON_GLA))*np.nan
    NMAD = np.zeros(len(bins_ON_GLA))*np.nan
    meas_area = np.zeros(len(bins_ON_GLA))*np.nan
    total_area = np.zeros(len(bins_ON_GLA))*np.nan

    # print 'Start bin loop'
    for i in np.arange(len(bins_ON_GLA)):
        # print 'loop nb: '+str(i)
        # t0=time.time()
        idx_alti = np.array(DEM>bins_ON_GLA[i]) & np.array(DEM<=bins_ON_GLA[i]+bin_alti)  & np.array(mask_ON_GLA ==1) #all the pixels ON GLA within the elev band
        Mean_alti[i] = bins_ON_GLA[i]+bin_alti/2.
        nb_pixel_band = len(DEM[idx_alti])
        # print time.time()-t0
        total_area[i] = nb_pixel_band*pixel_size*pixel_size
        dh_temp = dh_trend_ON_GLA[idx_alti].flatten()
        if len(dh_temp[~np.isnan(dh_temp)])>0:
            median_temp = np.nanmedian(dh_temp)
            Mediane_Dh_all[i] = median_temp
            MAD_temp = np.median(np.absolute(dh_temp[~np.isnan(dh_temp)]-median_temp))
            NMAD_temp = 1.4826*MAD_temp
            NMAD[i] = NMAD_temp
            # print time.time()-t0
            #dh_trend_ON_GLA[idx_alti & np.array(np.absolute(dh_trend_ON_GLA-median_temp)>3*NMAD_temp)] = np.nan
            dh_temp[np.absolute(dh_temp-median_temp)>3*NMAD_temp]=np.NaN
            #mean_band_check = np.nanmean(dh_trend_ON_GLA[idx_alti])
            mean_band = np.nanmean(dh_temp)
            Mean_Dh_3nmad[i] = mean_band
            Std_dh_3nmad[i] = np.nanstd(dh_temp)
            # print mean_band
            # print mean_band_check
            # print time.time()-t0
            #nb_pixel_meas_2 = len(dh_trend_ON_GLA[idx_alti & np.array(~np.isnan(dh_trend_ON_GLA)) ])
            nb_pixel_meas = len(dh_temp[~np.isnan(dh_temp)])
            # print nb_pixel_meas_2
            # print nb_pixel_meas
            meas_area[i] = nb_pixel_meas*pixel_size*pixel_size
            #dh_trend_ON_GLA[idx_alti & np.array(np.isnan(dh_trend_ON_GLA)) ] = mean_band_check #on remplit la bonde avec la valeure mediane
        # print time.time()-t0

    #dh_trend_ON_GLA = dh_trend_ON_GLA*ice_dens
    mean_all = np.nansum(Mean_Dh_3nmad*total_area*ice_dens)/np.nansum(total_area)
    #mean_all_2 = np.nanmean(dh_trend_ON_GLA[np.array(mask_ON_GLA >0)])
    #print 'Mass balance from mean: ',np.nanmean(dh_trend_ON_GLA[np.array(mask_ON_GLA >0)]), ' m w.e. a-1'
    print('Mass balance from mean:',mean_all,' m w.e. a-1')

    Mean_alti = np.array(Mean_alti)
    Mediane_Dh_all = np.array(Mediane_Dh_all)
    Mean_Dh_3nmad = np.array(Mean_Dh_3nmad)
    NMAD = np.array(NMAD)
    meas_area = np.array(meas_area)
    total_area = np.array(total_area)
    Std_dh_3nmad = np.array(Std_dh_3nmad)
    OUTPUT =  np.column_stack(( Mean_alti,Mediane_Dh_all,Mean_Dh_3nmad,NMAD,meas_area,total_area,Std_dh_3nmad))

    return mean_all, OUTPUT, dh_trend_ON_GLA

def bin_elev(dem_ref,mask,bin_type='fixed',bin_val=100.):

    dem_on_mask = np.copy(dem_ref)
    dem_on_mask [ np.invert(mask) ] = np.NaN

    min_elev = np.nanmin(dem_on_mask) - (np.nanmin(dem_on_mask) % bin_val)
    max_elev = np.nanmax(dem_on_mask) + 1

    if bin_type == 'fixed':
        bin_final = bin_val
    elif bin_type == 'percentage':
        bin_final = np.ceil(bin_val/100.*(max_elev-min_elev))
    else:
        sys.exit('Bin type not recognized.')

    bins_on_mask = np.arange(min_elev, max_elev, bin_final)

    a,b=np.shape(dem_on_mask)
    idx_bin = np.zeros((a,b,len(bins_on_mask)),dtype=np.bool)
    mean_bin = np.zeros(len(bins_on_mask))
    for i in np.arange(len(bins_on_mask)):
        idx_bin[:,:,i] = np.array(dem_on_mask>=bins_on_mask[i]) & np.array(dem_on_mask<(bins_on_mask[i]+bin_final))
        mean_bin[i] = bins_on_mask[i] + bin_final/2.

    return idx_bin, mean_bin

def mb_direct(dh_dt,mask,pixel_size):

    dh_dt_on_mask = np.copy(dh_dt)
    dh_dt_on_mask[mask != 1] = np.nan

    nb_pixel_total = len(mask[mask == 1])
    dh_all = np.nansum(dh_dt_on_mask)/nb_pixel_total
    total_area = nb_pixel_total * pixel_size * pixel_size

    return dh_all, total_area


def mb_hypso(dh_dt,dem_ref,mask,pixel_size,bin_type='fixed',bin_val=100.,filt_bin=None,interp_interbin=None,density='huss',voidfill=False):

    #reference DEM on mask only
    dem_on_mask = np.copy(dem_ref)
    dem_on_mask[mask != 1] = np.nan

    #dh on mask only
    dh_dt_on_mask = np.copy(dh_dt)
    dh_dt_on_mask[mask != 1] = np.nan

    #define elevation binning
    idx_bin, mean_bin = bin_elev(dem_ref,mask,bin_type,bin_val)
    nb_bin = len(mean_bin)

    #preallocate
    mean_dh_filtered = np.zeros(nb_bin)*np.nan
    std_dh_filtered = np.zeros(nb_bin)*np.nan
    nb_pixel_bin = np.zeros(nb_bin)*np.nan
    nb_pixel_novoid = np.zeros(nb_bin)*np.nan
    nb_pixel_filtered = np.zeros(nb_bin)*np.nan
    meas_area = np.zeros(nb_bin)*np.nan
    total_area = np.zeros(nb_bin)*np.nan

    if voidfill:
        dh_dt_voidfill = np.copy(dh_dt_on_mask)
    else:
        dh_dt_voidfill = None

    #loop for each elevation bin
    for i in np.arange(nb_bin):

        nb_pixel_bin[i] = len(dem_on_mask[idx_bin[:,:,i]])
        total_area[i] = nb_pixel_bin[i]*pixel_size*pixel_size
        dh_bin = dh_dt_on_mask[idx_bin[:,:,i]].flatten()

        nb_pixel_novoid[i] = len(dh_bin[~np.isnan(dh_bin)])

        if len(dh_bin[~np.isnan(dh_bin)])>0:

            if filt_bin == '3NMAD':
                median_temp = np.nanmedian(dh_bin)
                MAD_temp = np.nanmedian(np.absolute(dh_bin[~np.isnan(dh_bin)]-median_temp))
                NMAD_temp = 1.4826*MAD_temp
                dh_bin[np.absolute(dh_bin-median_temp)>3*NMAD_temp]=np.NaN
                if voidfill:
                    dh_dt_voidfill[idx_bin[:,:,i] & np.array(np.absolute(dh_dt_voidfill - median_temp) > 3 * NMAD_temp)] = np.nan

            mean_dh_filtered[i] = np.nanmean(dh_bin)
            std_dh_filtered[i] = np.nanstd(dh_bin)
            nb_pixel_filtered[i] = len(dh_bin[~np.isnan(dh_bin)])
            meas_area[i] = nb_pixel_filtered[i]*pixel_size*pixel_size


    #inter-bin interpolation
    if interp_interbin == 'mean':
        mean_all_bins = np.nanmean(mean_dh_filtered)
        mean_dh_filtered[np.isnan(mean_dh_filtered)] = mean_all_bins
    elif interp_interbin == 'pwlinear':
        idx_novoid = ~np.isnan(mean_dh_filtered)
        mean_bin_novoid = mean_bin[idx_novoid]
        dh_novoid = mean_dh_filtered[idx_novoid]
        if len(mean_bin_novoid) <= 1:
            print('Not enough non-void bins to interpolate: using mean to interpolate...')
            mean_all = np.nanmean(mean_dh_filtered)
            mean_dh_filtered[np.isnan(mean_dh_filtered)] = mean_all
        else:
            #interpolate linearly before first non void bin
            xp = mean_bin_novoid
            yp = dh_novoid
            x= mean_bin[~idx_novoid]

            mean_dh_filtered[~idx_novoid]=np.interp(x,xp,yp)
    elif interp_interbin == 'pwlinear2':
        idx_novoid = ~np.isnan(mean_dh_filtered)
        mean_bin_novoid = mean_bin[idx_novoid]
        dh_novoid = mean_dh_filtered[idx_novoid]
        if len(mean_bin_novoid) <= 1:
            print('Not enough non-void bins to interpolate: using mean to interpolate...')
            mean_all = np.nanmean(mean_dh_filtered)
            mean_dh_filtered[np.isnan(mean_dh_filtered)] = mean_all
        else:
            # interpolate linearly before first non void bin
            xp = mean_bin_novoid
            yp = dh_novoid
            x = mean_bin[~idx_novoid]
            mean_dh_filtered[~idx_novoid] = np.interp(x, xp, yp)

            #adjust values after and before
            if xp[0] != mean_bin[0]:
                a = (yp[1]-yp[0])/(xp[1]-xp[0])
                b = yp[1]- a * xp[1]

                idx_lowvoid = ~idx_novoid & np.array(mean_bin < xp[0])
                mean_dh_filtered[idx_lowvoid] = a * mean_bin[idx_lowvoid] + b

            if xp[-1] != mean_bin[-1]:
                a = (0-yp[-1])/(mean_bin[-1]-xp[-1])
                b = yp[-1]- a * xp[-1]

                idx_upvoid = ~idx_novoid & np.array(mean_bin > xp[-1])
                mean_dh_filtered[idx_upvoid] = a * mean_bin[idx_upvoid] + b
    elif interp_interbin == 'poly':
        idx_novoid = ~np.isnan(mean_dh_filtered)
        mean_bin_novoid = mean_bin[idx_novoid]
        dh_novoid = mean_dh_filtered[idx_novoid]
        if len(mean_bin_novoid) <= 1:
            print('Not enough non-void bins to interpolate: using mean to interpolate...')
            mean_all = np.nanmean(mean_dh_filtered)
            mean_dh_filtered[np.isnan(mean_dh_filtered)] = mean_all
        else:
            # interpolate linearly before first non void bin
            xp = mean_bin_novoid
            yp = dh_novoid
            x = mean_bin[~idx_novoid]

            P=np.polyfit(xp,yp,3)
            mean_dh_filtered[~idx_novoid] = np.polyval(P,x)
    elif interp_interbin == 'poly2':
        idx_novoid = ~np.isnan(mean_dh_filtered)
        mean_bin_novoid = mean_bin[idx_novoid]
        dh_novoid = mean_dh_filtered[idx_novoid]
        if len(mean_bin_novoid) <= 1:
            print('Not enough non-void bins to interpolate: using mean to interpolate...')
            mean_all = np.nanmean(mean_dh_filtered)
            mean_dh_filtered[np.isnan(mean_dh_filtered)] = mean_all
        else:
            # interpolate linearly before first non void bin
            xp = mean_bin_novoid
            xp.append(mean_bin[-1])
            yp = dh_novoid
            yp.append(0)
            x = mean_bin[~idx_novoid]

            P=np.polyfit(xp,yp,3)
            mean_dh_filtered[~idx_novoid] = np.polyval(P,x)

    #void filling
    if voidfill:
        for i in np.arange(nb_bin):
            dh_dt_voidfill[idx_bin[:,:,i] & np.array(np.isnan(dh_dt_voidfill))] = mean_dh_filtered[i]

    #glacier-wide volume change
    dh_all = np.nansum(mean_dh_filtered*total_area)/np.nansum(total_area)
    print('Volume change is:', dh_all, ' m a-1')

    #glacier-wide mass balance
    if density == 'huss':
        mb_all = dh_all * 0.85
    else:
        sys.exit('Density calculation not recognized')

    output =  np.column_stack(( mean_bin,mean_dh_filtered,meas_area,total_area,std_dh_filtered))

    return dh_all, mb_all , output, dh_dt_voidfill



#renewing stuff a bit in here:

def std_err_finite(std,Neff,neff):

    return std*np.sqrt(1/Neff*(Neff-neff)/Neff)

def std_err(std,Neff):

    return std*np.sqrt(1/Neff)

def linear_err(delta_x,std_acc_y):

    return delta_x**2/8*std_acc_y

def gauss(n=11,sigma=1):
    r = range(-int(n/2),int(n/2)+1)
    return [1 / (sigma * sqrt(2*pi)) * exp(-float(x)**2/(2*sigma**2)) for x in r]

def idx_near_val(array,v):

    return np.nanargmin(np.abs(array - v))

def interp_linear(xp,yp,errp,acc_y,loo=False):

    #interpolation 1d: nan are considered void, with possible leave-one-out (reinterpolate each value)

    #getting void index
    idx_void = np.isnan(yp)

    #preallocating arrays
    yp_out = np.copy(yp)
    errp_out = np.copy(errp)
    errlin_out = np.zeros(len(yp))*np.nan

    #don't really care about performance, let's do this one at a time
    for i in np.arange(len(xp)):

        x0=xp[i]
        tmp_xp = np.copy(xp)
        tmp_xp[idx_void] = np.nan
        if loo:
            tmp_xp[i] = np.nan #this is for leave-one out
        else:
            if not np.isnan(tmp_xp[i]):
                continue

        # find closest non void bin
        idx_1 = idx_near_val(tmp_xp, x0)
        tmp_xp[idx_1] = np.nan
        # second closest
        idx_2 = idx_near_val(tmp_xp, x0)

        #linear interpolation (or extrapolation)
        a = (xp[idx_2] - x0)/(xp[idx_2] - xp[idx_1])
        b = (x0 - xp[idx_1])/(xp[idx_2] - xp[idx_1])

        #propagating standard error
        y0_out = a * yp[idx_1] + b*yp[idx_2]
        err0_out = np.sqrt(a**2 * errp[idx_1]**2 + b**2 * errp[idx_2]**2)
        # err0_out = np.sqrt(errp[idx_1]**2 + errp[idx_2]**2)

        #estimating linear error
        delta_x = max(np.absolute(xp[idx_2] - x0),np.absolute(xp[idx_1]-x0))
        errlin0_out = linear_err(delta_x,acc_y)

        #appending
        yp_out[i] = y0_out
        errp_out[i] = err0_out
        errlin_out[i] = errlin0_out

    return yp_out, errp_out, errlin_out

def interp_lowess(xp,yp,errp,acc_y,rang,kernel='Exc'):


    yp_out, errp_out = lowess_homemade_kern(xp,yp,1/(errp**2),a1=rang/4.,kernel=kernel)

    idx_void = np.isnan(yp)
    errlin_out = np.zeros(len(yp))*np.nan

    for i in np.arange(len(xp)):

        x0=xp[i]
        tmp_xp = np.copy(xp)
        tmp_xp[idx_void] = np.nan
        if not np.isnan(tmp_xp[i]):
            continue

        # find closest non void bin
        idx_1 = idx_near_val(tmp_xp, x0)
        tmp_xp[idx_1] = np.nan

        delta_x = np.absolute(xp[idx_1] - x0)
        errlin0_out = linear_err(delta_x, acc_y)
        errlin_out[i] = np.sqrt(errlin0_out**2+errp[idx_1]**2)

    return yp_out, errp_out, errlin_out


def vol_hypso_linear(ddem, dem, mask, gsd, slope, neff_geo, neff_num, acc_dh, std_stable, rang, estim_std = None, bin_type='fixed',bin_val=50.,filt_bin='3NMAD',method='lowess'):

    #mask only valid pixels in the mask
    final_mask = np.logical_and(np.logical_and(np.isfinite(ddem), np.isfinite(dem)), mask)

    dem_on_mask = dem[final_mask]
    ddem_on_mask = ddem[final_mask]

    #for void-filled output
    ddem_out = np.copy(ddem)

    #binning
    min_elev = np.nanmin(dem[mask]) - (np.nanmin(dem[mask]) % bin_val)
    max_elev = np.nanmax(dem[mask]) + 1
    if bin_type == 'fixed':
        bin_final = bin_val
    elif bin_type == 'percentage':
        bin_final = np.ceil(bin_val / 100. * (max_elev - min_elev))
    else:
        sys.exit('Bin type not recognized.')
    bins_on_mask = np.arange(min_elev, max_elev, bin_final)

    #preallocating
    nb_bin = len(bins_on_mask)
    elev_bin = np.zeros(nb_bin) * np.nan
    nmad_bin = np.zeros(nb_bin) * np.nan
    mean_bin = np.zeros(nb_bin) * np.nan
    med_bin = np.zeros(nb_bin) * np.nan
    std_bin = np.zeros(nb_bin) * np.nan
    slope_bin = np.zeros(nb_bin) * np.nan
    area_tot_bin = np.zeros(nb_bin) * np.nan
    area_meas_bin = np.zeros(nb_bin) * np.nan
    std_err_bin = np.zeros(nb_bin) * np.nan
    std_fin_bin = np.zeros(nb_bin) * np.nan
    nonvoid_err_bin = np.zeros(nb_bin) * np.nan


    #do this one bin at a time to avoid filling in memory with a lot of masks
    for i in np.arange(nb_bin):

        idx_bin = np.array(dem_on_mask >= bins_on_mask[i]) & np.array(
            dem_on_mask < (bins_on_mask[i] + bin_final))
        idx_orig = np.array(dem >= bins_on_mask[i]) & np.array(
            dem < (bins_on_mask[i] + bin_final)) & mask
        area_tot_bin[i] = np.count_nonzero(idx_orig)*gsd**2
        area_meas_bin[i] = np.count_nonzero(idx_bin)*gsd**2
        elev_bin[i] = bins_on_mask[i] + bin_final / 2.
        dh_bin = ddem_on_mask[idx_bin]
        slope_bin[i] = np.nanmedian(slope[idx_orig])

        if len(dh_bin[~np.isnan(dh_bin)]) > 0:

            med_bin[i] = np.nanmedian(dh_bin)
            if filt_bin=='3NMAD':
                mad = np.nanmedian(np.absolute(dh_bin - med_bin[i]))
                nmad = 1.4826 * mad
                nmad_bin[i] = nmad
                idx_outlier = np.absolute(dh_bin - med_bin[i]) > 3*nmad
                nb_outlier = np.count_nonzero(idx_outlier)
                dh_bin[idx_outlier] = np.nan
                ddem_out[idx_orig & np.array(np.absolute(ddem_out - med_bin[i]) > 3 * nmad)] = np.nan
                area_meas_bin[i] -= nb_outlier*gsd**2
            std_bin[i] = np.nanstd(dh_bin)
            mean_bin[i] = np.nanmean(dh_bin)
            ddem_out[idx_orig & np.isnan(ddem_out)] = mean_bin[i]


    #first, get standard error for all non-void bins
    idx_nonvoid = area_meas_bin>0

    area_tot = np.sum(area_tot_bin)

    if estim_std is not None:
        std_bin = estim_std

    std_fin_bin[idx_nonvoid] = std_err_finite(std_bin[idx_nonvoid],neff_geo*area_tot_bin[idx_nonvoid]/area_tot,neff_geo*area_meas_bin[idx_nonvoid]/area_tot)
    std_err_bin[idx_nonvoid] = std_err(std_stable,neff_num*area_meas_bin[idx_nonvoid]/area_tot)
    nonvoid_err_bin[idx_nonvoid] = np.sqrt(std_fin_bin[idx_nonvoid]**2 + std_err_bin[idx_nonvoid]**2)

    if method == 'linear':
        #first, do a leave-one out linear interpolation to remove non-void bins with really low confidence
        loo_mean, loo_std_err, loo_lin_err = interp_linear(elev_bin,mean_bin,nonvoid_err_bin,acc_dh,loo=True)
        loo_full_err = np.sqrt(loo_std_err**2+loo_lin_err**2)

        idx_low_conf = nonvoid_err_bin>loo_full_err

        idx_final_void = np.logical_and(np.invert(idx_nonvoid),idx_low_conf)

        #then, interpolate for all of those bins
        mean_bin[idx_final_void]=np.nan
        nonvoid_err_bin[idx_final_void]=np.nan

        final_mean, final_std_err, final_lin_err = interp_linear(elev_bin,mean_bin,nonvoid_err_bin,acc_dh,loo=False)

        final_std_err[~idx_final_void] = 0

    elif method == 'lowess':

        final_mean, final_std_err, final_lin_err = interp_lowess(elev_bin,mean_bin,nonvoid_err_bin,acc_dh,rang)

        final_std_err[idx_nonvoid]=0

    else:
        print('Inter-bin interpolation method must be "linear" or "lowess"')
        sys.exit()

    final_std_err[np.isnan(final_std_err)] = 0
    final_lin_err[np.isnan(final_lin_err)] = 0
    interbin_err = np.sqrt(final_std_err**2+final_lin_err**2)
    intrabin_err = std_fin_bin
    intrabin_err[np.isnan(intrabin_err)] = 0
    final_mean[idx_nonvoid]=mean_bin[idx_nonvoid]

    tot_err = np.sqrt(interbin_err**2+intrabin_err**2)

    df = pd.DataFrame()
    df =df.assign(elev=elev_bin,mean_dh=mean_bin,std_dh=std_bin,slope=slope_bin,f_mean=final_mean,intra_err=intrabin_err,inter_err=interbin_err,area_tot=area_tot_bin,area_meas=area_meas_bin,tot_err=tot_err)

    for i in np.arange(nb_bin):
        idx_orig = np.array(dem >= bins_on_mask[i]) & np.array(
            dem < (bins_on_mask[i] + bin_final)) & mask
        if not idx_nonvoid[i]:
            ddem_out[idx_orig] = final_mean[i]

    return df, ddem_out


if __name__ == '__main__':

    fn_ddem = '/home/atom/ongoing/std_err/data_vhr/etienne_mb/dh_MB_FULL_10m.tif'
    fn_dem = '/home/atom/ongoing/std_err/data_vhr/etienne_mb/DEM_REF_MB_FULL_10m.tif'
    fn_maskvoid = '/home/atom/ongoing/std_err/data_vhr/etienne_mb/mask_abla_void2.tif'
    fn_mask = '/home/atom/ongoing/std_err/data_vhr/etienne_mb/mask_mdg.tif'

    from rastlib import read_nanarray

    ddem = read_nanarray(fn_ddem)
    dem = read_nanarray(fn_dem)
    maskvoid = (read_nanarray(fn_maskvoid) ==1)
    mask = (read_nanarray(fn_mask) ==1)

    ddem[np.absolute(ddem)>50] = np.nan
    ddem[~maskvoid] = np.nan

    gsd = 10.
    neff_geo=33
    neff_num=20000
    std_stable=1
    acc_dh= 1/300.
    bin_type='fixed'
    bin_val=50.
    filt_bin = '3NMAD'
    method='lowess'

    df, _ = vol_hypso_linear(ddem, dem, mask, gsd, neff_geo, neff_num, acc_dh, std_stable, var_range=100.)

    # from demlib import dem_contour_fl
    #
    # elev_contour = elev-25.
    # elev_contour = np.array(elev_contour,elev_contour[-1]+50.)
    # fn_shp_out = '/home/atom/ongoing/std_err/data_vhr/etienne_mb/elev_contour.shp'
    # dem_contour_fl(fn_dem,fn_shp_out,elev-25.)
    # from vectlib import clip_shp_to_shp
    #
    # fn_shp ='/home/atom/ongoing/Test_SGS/outlines/MdG.shp'
    # fn_clipped_out ='/home/atom/ongoing/std_err/data_vhr/etienne_mb/elev_contour_clipped.shp'
    # clip_shp_to_shp(fn_shp_out,fn_shp,fn_clipped_out)

    mask_void_2 = read_nanarray(fn_maskvoid)
    mask_void_2[~mask]=-9999

    from rastlib import write_nanarray
    write_nanarray('/home/atom/ongoing/std_err/figures/fig6/maskvoid.tif',fn_ddem,mask_void_2)

    out_dir='/home/atom/ongoing/std_err/figures/fig6'
    df.to_csv(os.path.join(out_dir,'df_interp_lowess.csv'))



