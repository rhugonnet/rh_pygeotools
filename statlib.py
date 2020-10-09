
import os, sys
import numpy as np
from scipy import stats, linalg
from math import ceil

def RegLinWeightedMat(x, y, w,conf_interv=0.99, conf_slope = 0.95):
    X = x*1.0
    Y = y *1.0
    W = w * 1.0

    Y[np.isnan(W) | np.isnan(X)] = np.nan #check for NaNs
    # X[np.isnan(W) | np.isnan(Y)] = np.nan #check for NaNs
    W[np.isnan(Y) | np.isnan(X)] = np.nan #check for NaNs

    sum_w = np.nansum(W, axis = 0)
    moy_X_w = np.nansum(X*W, axis = 0)/sum_w
    moy_Y_w = np.nansum(Y*W, axis = 0)/sum_w

    mat_cross_product = W*(X-moy_X_w)*(Y-moy_Y_w)
    sum_mat_cross_product = np.nansum(mat_cross_product, axis = 0)

    mat_X_squared = W*(X-moy_X_w)**2
    sum_mat_X_squared = np.nansum(mat_X_squared, axis = 0)

    beta1 = sum_mat_cross_product/sum_mat_X_squared
    beta0 = moy_Y_w - beta1*moy_X_w

    alpha_interv=1.-conf_interv
    alpha_slope = 1.-conf_slope

    Y_pred = beta1*X+beta0
    n = np.sum(~np.isnan(X),axis=0)
    SSX = sum_mat_X_squared
    SXY = np.sqrt(np.nansum(W*(Y-Y_pred)**2, axis = 0)/(n-2))
    SE_slope = SXY/np.sqrt(SSX)
    hi = 1./n+(X-moy_X_w)**2/SSX

    # quantile of student's t distribution for p=1-alpha/2
    q_interv = stats.t.ppf(1.-alpha_interv/2, n-2)
    q_slope = stats.t.ppf(1.-alpha_slope/2, n-2)

    # get the upper and lower CI:
    dy = q_interv*SXY*np.sqrt(hi)
    Yl = Y_pred-dy
    Yu = Y_pred+dy

    # calculate incert on slope
    incert_slope = q_slope*SE_slope

    return beta1, beta0, incert_slope, Yl, Yu


def kernel_exclass(xi,x0,a1,kappa=0.5):
    return np.exp(-(np.abs(xi - x0 )/a1) ** kappa)

def kernel_exp(xi,x0,a1):
    return np.exp(-np.abs(xi-x0)/a1)

def kernel_gaussian(xi,x0,a1):
    return np.exp(-((xi-x0)/a1)**2)

#TODO: kernel spherical?

def lowess_homemade_kern(x, y, w, a1,kernel='Exp'):
    """
    #inspired by: https://xavierbourretsicotte.github.io/loess.html

    homebaked lowess with variogram kernel + heteroscedasticity of observations with error

    :param x:
    :param y:
    :param w: heteroscedastic weights (inverse of variance)
    :param a1: range of the kernel (in variogram terms)
    :param kernel: kernel function
    :return:
    """

    n = len(x)
    yest = np.zeros(n)
    err_yest = np.zeros(n)

    if kernel == 'Gau':
        kernel_fun = kernel_gaussian
    elif kernel =='Exp':
        kernel_fun = kernel_exp
    elif kernel == 'Exc':
        kernel_fun = kernel_exclass
    else:
        print('Kernel not recognized.')
        sys.exit()

    # Initializing all weights from the bell shape kernel function
    W = np.array([kernel_fun(x,x[i],a1)*w for i in range(n)]).T

    # # Looping through all x-points
    # for i in range(n):
    #     weights = w[:, i]
    #     b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
    #     A = np.array([[np.sum(weights), np.sum(weights * x)],
    #                   [np.sum(weights * x), np.sum(weights * x * x)]])
    #     theta = linalg.solve(A, b)
    #     yest[i] = theta[0] + theta[1] * x[i]

    X = np.array([x for i in range(n)]).T
    Y = np.array([y for i in range(n)]).T

    beta1, beta0, _, Yl, Yu = RegLinWeightedMat(X,Y,W,conf_interv=0.68)

    for i in range(n):
        yest[i] = beta1[i]*x[i]+beta0[i]
        err_yest[i] = (Yu[i,i] - Yl[i,i])/2

    return yest, err_yest


def lowess_ag(x, y, f=2. / 3., iter=3):
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2

    return yest

x = np.linspace(0,1,100)
noise = np.random.normal(loc = 0, scale = .25, size = 100)
y = np.sin(x * 1.5 * np.pi )
y_noise = y + noise



if __name__ == '__main__':

    x=elev_bin
    y=mean_bin
    w=1/(nonvoid_err_bin**2)
    a1=200.