#!/home/fs71386/lkugler/miniconda3/bin/python
"""Convert University of Wyoming Database format to WRF's sounding format."""

import pandas as pd
import numpy as np
import csv, sys, os

def write_wrfprof(p_2m, z, pot_tmp, r, u, v, f_out):
    # surface measurements
    t_2m = pot_tmp[0]  # surface potential Temperature
    r_2m = r[0]  # surface vapor mixing ratio
    n_levels = z.shape[0]

    line1 = '{:9.2f} {:9.2f} {:10.2f}'.format(p_2m, t_2m, r_2m)
    wrfformat = '{:9.2f} {:9.2f} {:10.2f} {:10.2f} {:10.2f}'

    with open(f_out, 'w') as f:
        f.write(line1+' \n')
        for i in range(n_levels-1):  # ignore last level due to nan
            d = wrfformat.format(z[i], pot_tmp[i], r[i], u[i], v[i])
            #print(d)
            f.write(d+ '\n')
    print(f_out, 'saved.')


def read_uwyo(f_in):
    skiprows = 6  # make sure this is correct
    width = 7 # width of a data column
    columns = 11

    with open(f_in) as f:
        for i, line in enumerate(f):
            pass
    N_lines = i

    with open(f_in) as f:
        N_lines = N_lines-skiprows
        data = np.full((N_lines, 11), np.nan)
        j = 0

        for i, line in enumerate(f):
            if i>skiprows:
                elements = [line[i*width:(i+1)*width] for i in range(columns)]

                for i, element in enumerate(elements):
                    if element == 7*' ':
                        elements[i] = np.nan
                    else:
                        elements[i] = float(element.strip())

                data[j,:] = np.array(elements)
                j += 1
    print(data)

    p = data[:, 0]
    z = data[:, 1]
    T = data[:, 2]
    Td = data[:, 3]
    r = data[:, 5]
    r[np.isnan(r)] = 0.00
    dd = data[:, 6]
    ff = data[:, 7]
    pot_tmp = data[:, 8]

    return p, z, T, Td, r, dd, ff, pot_tmp



base = '/home/fs71386/lkugler/wrf_sounding/'
f_in = base+'data/uwyo/06610_20080730_12z.txt'

f_out = base+'data/wrf/06610_2008073012_uwyo.wrfprof'
f_out_ens = base+'data/wrf/ens/from_uwyo/06610_2008073012_uwyo.*.wrfprof'
os.makedirs(os.path.dirname(f_out_ens), exist_ok=True)

save_csv = False
write_ensemble = False
plot_ensembles = False
###########

p, z, T, Td, r, dd, ff, pot_tmp = read_uwyo(f_in)

#### convert to WRF variables
ff = ff * 1852/3600  # conversion from knots to meters/hour to meters/second
u = -ff * np.cos(dd/180*np.pi-np.pi/2)
v = -ff * np.sin(dd/180*np.pi+np.pi/2)


if save_csv:
    df = pd.DataFrame(data={'p': p, 'T': T, 'Qv': r})
    csvname = ''.join(f_in.split('.')[:-1])+'.csv'
    df.to_csv(csvname)
    print(csvname, 'saved.')


p_2m = p[0]  # surface pressure
write_wrfprof(p_2m, z, pot_tmp, r, u, v, f_out)


if write_ensemble:
    # add perturbations from f_pert and add it to original profile
    # write perturbed profile to another file

    if plot_ensembles:
        import matplotlib.pyplot as plt

    for i in range(41):
        f_pert = base+'data/LMU/pert/raso.raso.'+str(i).zfill(3)+'.csv'
        print(f_in, f_pert)
        f_out = f_out_ens.replace('*', str(i).zfill(3))
        df = pd.read_csv(f_pert)

        def interp(old_y):
            old_z = df['z']
            new_z = z
            return np.interp(new_z, old_z, old_y)

        per_t = interp(df['potT'].values)
        per_r = interp(df['Qv'].values)
        per_u = interp(df['u'].values)
        per_v = interp(df['v'].values)
        #print('perturbations:', per_t[0], per_r[0], per_u[0], per_v[0])

        pot_tmp0 = pot_tmp + per_t
        r0 = r + per_r
        u0 = u + per_u
        v0 = v + per_v

        if plot_ensembles:
            fig, ax = plt.subplots(1,2, figsize=(6,12), sharey=True)
            ax[0].plot(pot_tmp0, z)
            ax[1].plot(u0, z)
            ax[1].plot(v0, z)
            ax[1].barbs(np.zeros(len(z))[::2], z[::2], u0[::2], v0[::2])
            ax[1].set_title('u, v')
            ax[0].set_title('potential temperature')
            fig.savefig(f_out[:-7]+'png')
            print(f_out[:-7]+'png', 'saved.')
            plt.close()

        p_2m = p[0]  # surface pressure
        write_wrfprof(p_2m, z, pot_tmp0, r0, u0, v0, f_out)
