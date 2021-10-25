#!/home/fs71386/lkugler/miniconda3/envs/DART/bin/python
import pandas as pd
import numpy as np
import os, csv, glob, sys

import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

import plot_profile as prof

def convert_one(f_in, f_out=False, dir_out='./', save_csv=True, debug=False):
    """Convert LMU formatted profile to WRF formatted profile.
    Options: 
      save in csv format and plotting
    """
    # Read in LMU formatted profile
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
    ############################################

    p = data[:, 0]
    z = data[:, 1]
    T = data[:, 2]
    ff = data[:, 6]
    dd = data[:, 7]

    if True:
        # use RH column as input moisture profile
        rh = data[:, 4]/100
        svp = mpcalc.saturation_vapor_pressure(T*units.K)
        vp = rh * svp
        r = mpcalc.mixing_ratio(vp, p*units.millibar)
    else:
        # use the mixing ratio column as input moisture profile
        r = data[:, 5]/1000.

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
    # specify moisture profile (testing)
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

    ###############
    def plot(T, p, r, f_in, dir_out):
        T = T * units.K
        p = p * units.millibar

        vp = mpcalc.vapor_pressure(p, r)
        svp = mpcalc.saturation_vapor_pressure(T)
        Td = mpcalc.dewpoint_from_relative_humidity(T, vp/svp)
        Td[np.isnan(Td)] = -99.*units.degree_Celsius  # fill nan with very low temp

        f_in_basename = os.path.basename(f_in)
        prof.core(p, T, Td, u, v, 
                saveto=dir_out+'prof_'+f_in_basename+'.png', title=f_in_basename)

    plot(T, p, r, f_in, dir_out)

    ####################################
    # surface measurements
    sp = p[0]  # surface pressure
    t_2m = pot_tmp[0]  # surface potential Temperature
    r_2m = r[0].magnitude*1000.  # surface vapor mixing ratio
    n_levels = z.shape[0]

    line1 = '{:9.2f} {:9.2f} {:10.2f}'.format(sp, t_2m, r_2m)

    if save_csv and f_out:
        os.makedirs(dir_out+'/csv/', exist_ok=True)
        csvname = dir_out+'/'+'.'.join(f_out.split('.')[:-1])+'.csv'

        df = pd.DataFrame(data={'p': p, 'T': T, 'Qv': r})
        df.to_csv(csvname)
        print(csvname, 'saved.')

    wrfformat = '{:9.2f} {:9.2f} {:10.2f} {:10.2f} {:10.2f}'
    r *= 1000.

    if f_out:
        r = r.magnitude

        f_out = dir_out+f_out
        with open(f_out, 'w') as f:
            f.write(line1+' \n')
            for i in range(n_levels):
                d = wrfformat.format(z[i], pot_tmp[i], r[i], u[i], v[i])
                if debug: print(d)
                f.write(d+ '\n')
        print(f_out, 'saved.')

if __name__ == '__main__':
    # just plot
    if False:
        f_in = '/jetfs/home/lkugler/wrf_sounding/data/LMU/improved/raso.v2'
        convert_one(f_in, f_out=False)

    maindir = '/home/fs71386/lkugler/wrf_profiles/' #'/jetfs/home/lkugler/wrf_profiles/'

    infiles = glob.glob(maindir+'/data/LMU/improved_pert10/raso.fc.*')
    dir_out = maindir+'/data/wrf/ens/improved_pert10/'
    os.makedirs(dir_out, exist_ok=True)

    for f_in in infiles:
        fname_out = os.path.basename(f_in)+'.wrfprof'
        convert_one(f_in, fname_out, dir_out=dir_out)
        #sys.exit()
