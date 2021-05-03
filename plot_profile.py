#!/home/srvx11/lehre/users/a1254888/.conda/envs/WRF4/bin/python
import os, sys
import numpy as np
import datetime as dt
import pandas as pd
import xarray as xr
import wrf
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def core(p, T, Td, u, v, **kwargs):

    # Calculate the LCL
    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])

    #print('LCL p, t:', int(lcl_pressure), int(lcl_temperature))

    # Calculate the parcel profile.
    parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')


    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(8, 8))
    skew = SkewT(fig, rotation=45)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    #skew.plot(p, T, 'k-')
    skew.plot(p, T, 'r.-', ms=5, lw=2, label='mean T')
    skew.plot(p, Td, 'g.-', ms=5, lw=2, label='mean Td')
    skew.plot_barbs(p, u, v)
    skew.ax.set_ylim(1000, 180)
    skew.ax.set_xlim(-20, 40)

    # Plot LCL temperature as black dot
    skew.plot(lcl_pressure, lcl_temperature, 'k.', markerfacecolor='black')

    # Plot the parcel profile as a black line
    skew.plot(p, parcel_prof, 'k', linewidth=2)

    # Shade areas of CAPE and CIN
    skew.shade_cin(p, T, parcel_prof)
    skew.shade_cape(p, T, parcel_prof)

    # Plot a zero degree isotherm
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

    # Add the relevant special lines
    skew.plot_dry_adiabats(lw=.5)
    skew.plot_moist_adiabats(lw=.5)
    skew.plot_mixing_lines(lw=.5)

    # Show the plot
    #plt.show()
    #skew.ax.set_title(time_str)
    plt.legend(loc='lower left')
    plt.title(kwargs.get('title'))
    fname = kwargs.get('saveto', 'profile.png')
    fig.savefig(fname)
    print(fname, 'saved.')
    plt.close()



def core_ens(p, T, Td, u, v, p2, T2, Td2, **kwargs):
    #from IPython import embed; embed()
    T = T.to(units.K)
    Td = Td.to(units.K)

    Tmean = T.mean(axis=0)
    Tdmean = Td.mean(axis=0)
    pmean = p.mean(axis=0)
    umean, vmean = u.mean(axis=0), v.mean(axis=0)
    Tstd = np.std(T.data, axis=0) * units.K
    Tdstd = np.std(Td.data, axis=0) * units.K


    # Calculate the parcel profile.
    parcel_prof = mpcalc.parcel_profile(pmean, Tmean[0], Tdmean[0]).to('degC')


    # Create a new figure. The dimensions here give a good aspect ratio
    fig = plt.figure(figsize=(8, 8))
    skew = SkewT(fig, rotation=45)

    # Plot a zero degree isotherm
    skew.ax.axvline(0, color='k', linestyle='--', linewidth=1)

    # Add the relevant special lines
    skew.plot_dry_adiabats(lw=.5)
    skew.plot_moist_adiabats(lw=.5)
    skew.plot_mixing_lines(lw=.5)

    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    skew.shade_area(pmean, Tmean-Tstd, Tmean+Tstd, color='r', label='std$_{ens}$ mean$_{dom}$ T$_{fc}$')
    skew.shade_area(pmean, Tdmean-Tdstd, Tdmean+Tdstd, color='g', label='std$_{ens}$ mean$_{dom}$ Td$_{fc}$')
    #skew.plot(p, Tmean+np.std(T), '-', color='grey', lw=1, label='p.999(T)')
    skew.plot(pmean, Tmean, 'r:', ms=3, lw=1, label='mean$_{ens, dom}$ T$_{fc}$')
    skew.plot(pmean, Tdmean, 'g:', ms=3, lw=1, label='mean$_{ens, dom}$ Td$_{fc}$')
    skew.plot_barbs(pmean, umean, vmean)

    # Plot the parcel profile as a black line
    skew.plot(pmean, parcel_prof, 'k', linewidth=.5)

    # Shade areas of CAPE and CIN
    skew.shade_cin(pmean, Tmean, parcel_prof, alpha=0.2)
    skew.shade_cape(pmean, Tmean, parcel_prof, alpha=0.2)


    # nature
    skew.plot(p2, T2, 'r.-', ms=5, lw=2, label='mean$_{ens, dom}$ T$_{nature}$')
    skew.plot(p2, Td2, 'g.-', ms=5, lw=2, label='mean$_{ens, dom}$ Td$_{nature}$')


    skew.ax.set_ylim(1000, 180)
    skew.ax.set_xlim(-20, 40)

    plt.legend(loc='lower left')
    plt.title(kwargs.get('title'))
    fname = kwargs.get('saveto', 'profile.png')
    fig.savefig(fname)
    print(fname, 'saved.')
    plt.close()

