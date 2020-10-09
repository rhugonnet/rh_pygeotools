import numpy as np
from scipy import integrate
from rastlib import proximity_rast_ds, create_mem_raster_on_ref, read_nanarray, pixel_size

def neff_sum_circular(area,crange1,psill1,crange2,psill2):

    #short range variogram
    c1 = psill1 # partial sill
    a1 = crange1  # short correlation range

    #long range variogram
    c1_2 = psill2
    a1_2 = crange2 # long correlation range

    h_equiv = np.sqrt(area / np.pi)

    #hypothesis of a circular shape to integrate variogram model
    if h_equiv > a1_2:
        std_err = np.sqrt(c1 * a1 ** 2 / (5 * h_equiv ** 2) + c1_2 * a1_2 ** 2 / (5 * h_equiv ** 2))
    elif (h_equiv < a1_2) and (h_equiv > a1):
        std_err = np.sqrt(c1 * a1 ** 2 / (5 * h_equiv ** 2) + c1_2 * (1-h_equiv / a1_2+1 / 5 * (h_equiv / a1_2) ** 3))
    else:
        std_err = np.sqrt(c1 * (1-h_equiv / a1+1 / 5 * (h_equiv / a1) ** 3) + c1_2 * (1-h_equiv / a1_2+1 / 5 * (h_equiv / a1_2) ** 3))

    return std_err

def neff_rect_sum(area,width,crange1,psill1,crange2,psill2):

    def hcov_sum(h,crange1=crange1,psill1=psill1,crange2=crange2,psill2=psill2):

        return h*(cov(h,crange1,model='Sph',psill=psill1)+cov(h,crange2,model='Sph',psill=psill2))

    width = min(width,area/width)

    full_int = integrate_fun(hcov_sum,0,min(width,crange2))[0]
    if psill2 > width:
        bin_int = np.linspace(width,min(area/width,psill2),100)
        for i in range(len(bin_int)-1):
            low = bin_int[i]
            upp = bin_int[i+1]
            mid = bin_int[i] + (bin_int[i+1]- bin_int[i])/2
            piec_int = integrate_fun(hcov_sum, low, upp)[0]
            full_int += piec_int * 2/np.pi*np.arctan(width/(2*mid))

    std_err = np.sqrt(2*np.pi*full_int / area)

    return std_err


def neff(mask,fn_ref,crange,model='Sph',psill=1.,kappa=1/2,nugget=0):

    def hcov(h,crange=crange,model=model,psill=psill,kappa=kappa,nugget=nugget):

        return h*cov(h,crange,model=model,psill=psill,kappa=kappa,nugget=nugget)

    # fn_ref = '/home/atom/ongoing/std_err/data_vhr/etienne_mb/mask_mdg.tif'
    # mask = (read_nanarray(fn_ref) == 1)

    #get inside proximity of mask
    ds = create_mem_raster_on_ref(fn_ref)
    ds.GetRasterBand(1).WriteArray(mask)
    ds_prox = proximity_rast_ds(ds,val=0)
    prox = ds_prox.GetRasterBand(1).ReadAsArray()

    #pole of innacessibility
    poi_tup = np.where(prox==np.nanmax(prox))
    poi_x = poi_tup[0][0]
    poi_y = poi_tup[1][0]

    #get proximity from the pole
    mask_poi = np.zeros(np.shape(mask),dtype=bool)
    mask_poi[poi_x,poi_y] = 1
    ds.GetRasterBand(1).WriteArray(mask_poi)
    ds_prox_poi = proximity_rast_ds(ds,val=1)
    prox_poi = ds_prox_poi.GetRasterBand(1).ReadAsArray()

    #maximum distance in the mast from the poi
    max_prox = np.nanmax(prox_poi[mask])
    res = pixel_size(fn_ref)
    bin_size = 3 * res

    #bin with circular area
    bins_area = np.arange(0,max_prox,bin_size)
    tot_area = np.zeros(len(bins_area))
    obs_area = np.zeros(len(bins_area))

    for i in np.arange(len(bins_area)):

        tot_area[i] = np.count_nonzero((prox_poi >= bins_area[i]) & (prox_poi < bins_area[i] + bin_size))
        obs_area[i] = np.count_nonzero((prox_poi[mask] >= bins_area[i]) & (prox_poi[mask] < bins_area[i] + bin_size))

    frac_area = obs_area/tot_area

    #integrate piecewise:
    full_int = 0
    for i in np.arange(len(bins_area)):
        low= bins_area[i]
        upp= bins_area[i] + bin_size
        piec_int = integrate_fun(hcov,low,upp)[0]
        print(piec_int)
        full_int += piec_int*frac_area[i]

    area_tot = np.count_nonzero(mask)*res**2

    gamma_unity = full_int/area_tot

    return 1./(1-gamma_unity)


def integrate_fun(fun,low_limit,up_limit):

    return integrate.quad(fun,low_limit,up_limit)

def cov(h,crange,model='Sph',psill=1.,kappa=1/2,nugget=0):

    return (nugget + psill) - vgm(h,crange,model=model,psill=psill,kappa=kappa)

def vgm(h,crange,model='Sph',psill=1.,kappa=1/2,nugget=0):

    c0 = nugget #nugget
    c1 = psill #partial sill
    a1 = crange #correlation range
    s = kappa #smoothness parameter for Matern class

    if model == 'Sph':  # spherical model
        if h < a1:
            vgm = c0 + c1 * (3 / 2 * h / a1-1 / 2 * (h / a1) ** 3)
        else:
            vgm = c0 + c1
    elif model == 'Exp':  # exponential model
        vgm = c0 + c1 * (1-np.exp(-h / a1))
    elif model == 'Gau':  # gaussian model
        vgm = c0 + c1 * (1-np.exp(- (h / a1) ** 2))
    elif model == 'Exc':  # stable exponential model
        vgm = c0 + c1 * (1-np.exp(-(h/ a1)**s))

    return vgm