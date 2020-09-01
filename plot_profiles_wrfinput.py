#!/home/srvx11/lehre/users/a1254888/.conda/envs/WRF4/bin/python
import os, sys, glob
import numpy as np
import xarray as xr
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import wrf
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

basetemp = 300
nx = 200
areadims = ['south_north', 'west_east']

def wrf_get_date_str(ds):
    try:
        return [str(a.dt.strftime('%Y-%m-%d %H:%M').values) for a in ds.XTIME]
    except Exception as e:
        return ['',]


def wrf_get_datetime(ds):
    strings = wrf_get_date_str(ds)
    return pd.to_datetime(strings)

def listdir(p):
    l = os.listdir(p)
    return [a for a in l if not a.startswith('.')]


def plot_profile(f, fig, skew):

    ds = xr.open_dataset(f)
    ds = ds.isel(Time=0)

    Theta_m = ds.T.mean(areadims) + basetemp
    pm = (ds.P + ds.PB).mean(areadims)
    Qvm = (ds.QVAPOR).mean(areadims)

    U = wrf.destagger(ds.U, 2, meta=True)
    V = wrf.destagger(ds.V, 1, meta=True)

    um = U.mean(areadims)
    vm = V.mean(areadims)

    Tm = Theta_m*(pm/1e5)**(2/7)

    p = pm.values/100. * units.millibar
    Tm = Tm.values * units.K
    Qv = Qvm.values

    vp = mpcalc.vapor_pressure(p, Qv)
    svp = mpcalc.saturation_vapor_pressure(Tm)
    #svp[svp<1e-16*units.millibar] = 1e-16 * units.millibar

    Td = mpcalc.dewpoint_from_relative_humidity(Tm, vp/svp)
    Td[np.isnan(Td)] = -99.*units.degree_Celsius  # fill nan with very low temp
    # print('Dewpoints: ', Td)
    # print('Temperatures: ', Tm)

    u = um.values
    v = vm.values

    skew.plot(p, Tm, 'r.-', ms=1, lw=.5, label='mean T')
    skew.plot(p, Td, 'g.-', ms=5, lw=2 , label='mean Td')
    #skew.plot_barbs(p, u, v)
    return Tm, Td, p

####################
if __name__ == '__main__':

    expname = 'exp_v1.10_LMU+shear_filter'
    n_ens = 20
    filepattern = '/home/fs71386/lkugler/data/sim_archive/'+expname+'/wrfinput/*/wrfinput_d01'
    wrfinput_files = glob.glob(filepattern)

    fig = plt.figure(figsize=(8, 8))
    skew = SkewT(fig, rotation=45)

    # Add the relevant special lines
    skew.plot_dry_adiabats(lw=.5)
    skew.plot_moist_adiabats(lw=.5)
    skew.plot_mixing_lines(lw=.5)

    Tdata = np.zeros((n_ens, 50))+np.nan
    Tddata = np.zeros((n_ens, 50))+np.nan
    pdata = np.zeros((n_ens, 50))+np.nan

    for i, f in enumerate(wrfinput_files):
        print('ens', f)
        T, Td, p = plot_profile(f, fig, skew)

        Tdata[i, :] = T.astype(np.float32)
        pdata[i, :] = p.astype(np.float32)
        Tddata[i, :] = Td.astype(np.float32)
        #Tdata[1, i, :] = Td

    skew.ax.set_ylim(1000, 180)
    skew.ax.set_xlim(-18, 29)
    skew.ax.set_ylabel('pressure [hPa]')
    skew.ax.set_xlabel('temperature [$^\circ$C]')


    handles = [mpl.lines.Line2D([], [], color='r', marker='.', ms=1, lw=.5),
               mpl.lines.Line2D([], [], color='g', marker='.', ms=5, lw=2)]
    labels = ['Temperature', 'Dewpoint', ]
    plt.legend(handles, labels, loc='lower left')
    fig.savefig('wrfinput_'+expname+'.png', dpi=300)


    fig, ax = plt.subplots()

    Tmean = np.nanmean(Tdata, axis=0)
    pmean = np.nanmean(pdata, axis=0)
    #print(Tdata, Tmean)
    tpert = Tdata - Tmean
    from sklearn.decomposition import PCA
    tpert_pca = PCA(n_components=4).fit_transform(tpert.T)

    for i in range(4):
        ax.plot(tpert_pca[:,i], pmean)

    ax.invert_yaxis()
    ax.set_yscale('log')
    ax.set_ylabel('pressure [hPa]')
    ax.set_xlabel('temperature perturbation [K]')
    h = ax.plot([], [], 'k-')
    #plt.legend(h)
    fig.savefig('tpert.png', dpi=300)




    #######################
    fig = plt.figure(figsize=(8, 8))
    skew = SkewT(fig, rotation=45)

    # Add the relevant special lines
    skew.plot_dry_adiabats(lw=.5)
    skew.plot_moist_adiabats(lw=.5)
    skew.plot_mixing_lines(lw=.5)

    Tdmean = np.nanmean(Tddata, axis=0)

    skew.plot(p, Tmean-273.15, 'r.-', ms=5, lw=2, label='mean T')
    skew.plot(p, Tdmean, 'g.-', ms=5, lw=2 , label='mean Td')
    skew.ax.set_ylim(1000, 180)
    skew.ax.set_xlim(-18, 29)
    skew.ax.set_ylabel('pressure [hPa]')
    skew.ax.set_xlabel('temperature [$^\circ$C]')


    handles = [mpl.lines.Line2D([], [], color='r', marker='.', ms=5, lw=2),
               mpl.lines.Line2D([], [], color='g', marker='.', ms=5, lw=2)]
    labels = ['Temperature', 'Dewpoint', ]
    plt.legend(handles, labels, loc='lower left')
    fig.savefig('wrfinput_mean_'+expname+'.png', dpi=300)
