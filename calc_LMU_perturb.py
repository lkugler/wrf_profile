#!/home/fs71386/lkugler/miniconda3/bin/python
"""Calculate differences btw different files in LMU format and plot them"""

import pandas as pd
import numpy as np
import os, csv, glob, sys
import matplotlib.pyplot as plt

def read_LMU(f_in):
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
            #print(line)
            if i==skiprows:
                header=line
            if i>skiprows:
                values = line.split()
                #print(values)
                values = [float(v) for v in values]
                #print(values)
                data[j,:] = np.array(values)
                j += 1

    p = data[:, 0]
    z = data[:, 1]
    T = data[:, 2]
    r = data[:, 5]
    ff = data[:, 6]
    dd = data[:, 7]

    # constants (amsglossary)
    Rd = 287.05
    Rv = 461.51
    cpd = 1005.7  # specific heat dry air
    cpv = 1847.  # specific heat vapor
    kappa_moist = Rd/cpd*(1+r*Rv/Rd)/(1+r*cpv/cpd)

    # Temperature to potential temperature
    T = T*(1000/p)**kappa_moist
    u = -ff * np.cos(dd/180*np.pi-np.pi/2)
    v = -ff * np.sin(dd/180*np.pi+np.pi/2)

    return p, z, T, r, u, v


def calc_pert(f_in_ref, f_in, csvname):

    p, z, T, r, u, v = read_LMU(f_in_ref)
    p1, z1, T1, r1, u1, v1 = read_LMU(f_in)
    print(f_in_ref, f_in)

    # fig, ax = plt.subplots()
    # ax.plot(u, z)
    # ax.plot(u1, z)
    # fig.savefig('test.png')
    # sys.exit()
    # for i in range(len(u)):
    #     print(u[i], u1[i], '#', v[i], v1[i])
    #sys.exit()
    #p = p1-p
    T = T1-T
    #z = z1-z
    r = r1-r
    u = u1-u
    v = v1-v
    #print(z[u>5], u[u>5], v[u>5])

    df = pd.DataFrame(data={'z': z, 'potT': T, 'Qv': r, 'u': u, 'v': v})
    df.to_csv(csvname)
    print(csvname, 'saved.')

    fig, ax = plt.subplots(1,2, sharey=True)
    ax[0].plot(df['potT'], df['z'])
    ax[0].set_title('potT perturbations')
    ax[1].plot(df['u'], df['z'])
    ax[1].plot(df['v'], df['z'])
    ax[1].set_title('u, v perturbations')
    fig.savefig(csvname[:-4]+'.png')
    print(csvname[:-4]+'.png', 'saved.')
    plt.close()


if __name__ == '__main__':

    dir_out = '/home/fs71386/lkugler/wrf_sounding/data/LMU/pert/'
    f_in_ref = '/home/fs71386/lkugler/wrf_sounding/data/LMU/raso.raso'

    infiles = glob.glob('/home/fs71386/lkugler/wrf_sounding/data/LMU/raso.raso.0*')
    for f_in in infiles:

        f_out = dir_out+os.path.basename(f_in)+'.csv'
        calc_pert(f_in_ref, f_in, f_out)
