#!/home/fs71386/lkugler/miniconda3/envs/DART/bin/python
"""Create modified profiles"""

import pandas as pd
import numpy as np
import os, csv, glob, sys

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

def convert_one(f_in, f_out, save_csv=True, debug=False):
    skiprows = 3  # make sure this is correct

    with open(f_in) as f:
	    for i, line in enumerate(f):
		    pass
    N_lines = i

    with open(f_in) as f:
	    N_lines = N_lines-skiprows
	    data = np.full((N_lines, 8), np.nan)
	    j = 0

	    for i, line in enumerate(f):
		    if i==skiprows:
			    header=line
		    if i>skiprows:
			    values = line.split()
			    values = [float(v) for v in values]
			    #print(values)
			    data[j,:] = np.array(values)
			    j += 1

    if debug: print('input header:', header)

    p = data[:, 0]
    z = data[:, 1]
    T = data[:, 2]
    r = data[:, 5]/1000.
    ff = data[:, 6]
    dd = data[:, 7]

    #### convert to WRF variables

    # constants (amsglossary)
    Rd = 287.05
    Rv = 461.51
    cpd = 1005.7  # specific heat dry air
    cpv = 1847.  # specific heat vapor
    kappa_moist = Rd/cpd*(1+r*Rv/Rd)/(1+r*cpv/cpd)
    print('kappa_moist mean, std:', kappa_moist.mean(), kappa_moist.std())

    ####################################
    # Wind modifications
    ff = np.concatenate([np.linspace(2, 30, 36),
                         np.linspace(30, 0, len(p)-36) ])
    dd = np.concatenate([np.linspace(90, 210, 10),
                   np.linspace(213, 270, 20),
                   np.linspace(273, 360, len(p)-30 )])

    ####################################

    # Temperature to potential temperature
    pot_tmp = T*(1000/p)**kappa_moist

    u = -ff * np.cos(dd/180*np.pi-np.pi/2)
    v = -ff * np.sin(dd/180*np.pi+np.pi/2)

    ####################################
    # modifications
    #vp = mpcalc.vapor_pressure(p*units.millibar, r)
    if False:
        svp = mpcalc.saturation_vapor_pressure(T*units.K)
        #Td = mpcalc.dewpoint_from_relative_humidity(T*units.K, vp/svp)
        rh = 0.9
        vp = rh*svp
        r1 = mpcalc.mixing_ratio(vp, p*units.millibar)

        fig, ax = plt.subplots()
        ax.plot(r, p, 'r-')
        r[(r<r1)&(T>263)] = r1[(r<r1)&(T>263)]  # allow r to be greater than 90%

        ax.plot(r, p, 'g-')
        ax.invert_yaxis()
        fig.savefig('test.png')
    # vp = mpcalc.vapor_pressure(p*units.millibar, r)
    # Td = mpcalc.dewpoint_from_relative_humidity(T*units.K, vp/svp)
    # ax.plot(Td, p, 'g-')
    # ax.invert_yaxis()
    # fig.savefig('test.png')
    #
    # fig, ax = plt.subplots()
    # ax.plot(dd, p)
    # ax.invert_yaxis()
    # fig.savefig('test2.png')
    ####################################
    # surface measurements
    sp = p[0]  # surface pressure
    t_2m = pot_tmp[0]  # surface potential Temperature
    r_2m = r[0]*1000.  # surface vapor mixing ratio
    n_levels = z.shape[0]

    line1 = '{:9.2f} {:9.2f} {:10.2f}'.format(sp, t_2m, r_2m)

    if save_csv:
        df = pd.DataFrame(data={'p': p, 'T': T, 'Qv': r})
        csvname = '.'.join(f_out.split('.')[:-1])+'.csv'
        df.to_csv(csvname)
        print(csvname, 'saved.')
    # print(df.columns)
    # d = df[['Z [m],']].to_dict()
    wrfformat = '{:9.2f} {:9.2f} {:10.2f} {:10.2f} {:10.2f}'
    r *= 1000.

    with open(f_out, 'w') as f:
        f.write(line1+' \n')
        for i in range(n_levels):
            d = wrfformat.format(z[i], pot_tmp[i], r[i], u[i], v[i])
            if debug: print(d)
            f.write(d+ '\n')
    print(f_out, 'saved.')

if __name__ == '__main__':

    infiles = glob.glob('/home/fs71386/lkugler/wrf_sounding/data/LMU/raso.*')
    dir_out = '/home/fs71386/lkugler/wrf_sounding/data/wrf/ens/LMU+shear/'
    os.makedirs(dir_out, exist_ok=True)

    for f_in in infiles:
        f_out = dir_out+os.path.basename(f_in)+'.wrfprof'
        convert_one(f_in, f_out)
        #sys.exit()
