#!/home/fs71386/lkugler/miniconda3/bin/python
"""Convert `LMU` format to WRF's sounding format."""

import pandas as pd
import numpy as np
import os, csv, glob


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

    # Temperature to potential temperature
    pot_tmp = T*(1000/p)**kappa_moist

    u = -ff * np.cos(dd/180*np.pi-np.pi/2)
    v = -ff * np.sin(dd/180*np.pi+np.pi/2)

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
    for f_in in infiles:

        dir_out = '/home/fs71386/lkugler/wrf_sounding/data/wrf/ens/from_LMU/'
        f_out = dir_out+os.path.basename(f_in)+'.wrfprof'
        convert_one(f_in, f_out)
