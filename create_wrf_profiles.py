import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import metpy.calc as mpcalc
#from metpy.plots import Hodograph, SkewT
from metpy.units import units
from create_ens_profile_from_wrfinput import perturb

"""
Authors
J. Schröttle/ K. Bachmann: initial authors
L. Kugler: 
- New structure
- fixed bug in arctan
- only wind speed was perturbed, not direction
"""

# constants (amsglossary)
Rd = 287.05
Rv = 461.51
cpd = 1005.7  # specific heat dry air
cpv = 1847.  # specific heat vapor


def load_reference_profile(f_in = '/jetfs/home/lkugler/wrf_profiles/raso.v2'):
    a=np.loadtxt(f_in, skiprows=4)

    #retrieve profiles
    z=a[:,1]  # height / m
    p=a[:,0]  # pressure / hPa
    t=a[:,2]  # temperature / K
    # td=a[:,3] # empty field
    rh=a[:,4] # relative humidity / %
    # r=a[:,5]  # water vapor mixing ratio / (g/kg)
    ws=a[:,6] # wind speed / (m/s)
    wd=a[:,7] # wind direction / (deg)

    # nz = z.shape[0]  # number of layers
    return z, p, t, rh, ws, wd

def convert_natural_to_cartesian(ws, wd):
    u, v = mpcalc.wind_components(ws*units.mps, wd*units.degrees)
    return u.to('meter_per_second').m, v.to('meter_per_second').m

def convert_cartesian_to_natural(u, v):
    ws = np.sqrt(u**2 + v**2)
    return ws, mpcalc.wind_direction(u*units.mps,v*units.mps).to('degrees').m


def perturb_profile(z, t, rh, ws, wd):
    """Get perturbed vectors of temperature, relative humidity, ...
    
    Args:
        np.array
    Returns:
        np.array
    """
    step=20
    nz = len(z)
    nz_step = int(nz/step)
    print('number of random numbers:', nz_step)
    
    t += perturb(0., .25, nz_step, z)
    rh += rh*perturb(0, .02, nz_step, z)

    u, v = convert_natural_to_cartesian(ws, wd)
    u += perturb(0., .25, nz_step, z)
    v += perturb(0., .25, nz_step, z)
    # ws, wd = convert_cartesian_to_natural(u, v)

    return t, rh, u, v

def specify_wind_profile(p):
    # Wind modifications
    ws = np.concatenate([np.linspace(2, 30, 36),
                         np.linspace(30, 0, len(p)-36) ])
    wd = np.concatenate([np.linspace(90, 210, 10),
                   np.linspace(213, 270, 20),
                   np.linspace(273, 360, len(p)-30 )])
    return ws, wd

def calc_mixing_ratio(p, t, rh):

    rh = rh/100
    svp = mpcalc.saturation_vapor_pressure(t*units.K)
    vp = rh * svp
    r = mpcalc.mixing_ratio(vp, p*units.millibar)
    return r


def write_WRF_format(z, p_surface, pot_tmp, r, u, v, f_out='./test.txt'):
    # from IPython import embed; embed()
    
    if hasattr(r, 'magnitude'):
        r = r.to('g/kg').magnitude  
    if hasattr(z, 'magnitude'):
        z = z.magnitude
    if hasattr(p_surface, 'magnitude'):
        p_surface = p_surface.magnitude
    if hasattr(pot_tmp, 'magnitude'):
        pot_tmp = pot_tmp.magnitude
    if hasattr(u, 'magnitude'):
        u = u.magnitude
    if hasattr(v, 'magnitude'):
        v = v.magnitude
        
    
    # surface measurements
    sp = float(p_surface)  # surface pressure
    t_2m = float(pot_tmp[0])  # surface potential Temperature
    r_2m = float(r[0])  # surface vapor mixing ratio
    n_levels = z.shape[0]
    

    line1 = '{:9.2f} {:9.2f} {:10.2f}'.format(sp, t_2m, r_2m)
    wrfformat = '{:9.2f} {:9.2f} {:10.2f} {:10.2f} {:10.2f}'
    
    os.makedirs(os.path.dirname(f_out), exist_ok=True)
    with open(f_out, 'w') as f:
        f.write(line1+' \n')
        for i in range(n_levels):

            d = wrfformat.format(float(z[i]), float(pot_tmp[i]), float(r[i]), float(u[i]), float(v[i]))
            f.write(d+ '\n')
    print(f_out, 'saved.')



def profile(p, T, r, u, v,  f_out='./test.png'):
    T = T * units.K
    p = p * units.millibar

    vp = mpcalc.vapor_pressure(p, r)
    svp = mpcalc.saturation_vapor_pressure(T)
    Td = mpcalc.dewpoint_from_relative_humidity(T, vp/svp)
    Td[np.isnan(Td)] = -99.*units.degree_Celsius  # fill nan with very low temp
    import osselyze.plot_profile as pprof
    pprof.core(p, T, Td, u, v, figsize=(5,5.5), dpi=100, saveto=f_out)


def hodograph(z, u, v, f_out='./test.png'):

    def colored_line(x, y, z, ax, norm, lw=1, cmap='viridis'):
        """Line colored depending on z"""
        # add colorbar with `fig.colorbar(line, ax=ax)`
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        # Set the values used for colormapping
        lc.set_array(z)
        lc.set_linewidth(lw)
        line = ax.add_collection(lc)
        return line

    import ultraplot as plot
    fig, ax = plot.subplots(figsize=(3,3), journal='ams1')
    ax.set_aspect(1)
    for radius in [5, 10, 15, 20, 25, 30]:
        ax.add_patch(plt.Circle((0, 0), radius, color='grey', fill=False))

    norm = plt.Normalize(1, 12)
    line = colored_line(u, v, z/1000, ax, norm, lw=3)

    cb = ax.colorbar(line, loc='r', label='height [km]', values=[1, 5, 10, 15, 20])

    mmin, mmax = min(u.min(), v.min())-1, max(u.max(), v.max())+1
    ax.plot([0, 0], [mmin, mmax], ls='--', lw=1, color='grey')
    ax.plot([mmin, mmax], [0, 0], ls='--', lw=1, color='grey') 
    ax.set_xlim(mmin, mmax) #u.min()-1, u.max()+1)
    ax.set_ylim(mmin, mmax) #v.min()-1, v.max()+1)
    ax.set_xlabel('u [m/s]')
    ax.set_ylabel('v [m/s]')

    fig.savefig(f_out, transparent=False, dpi=200)
    print(f_out, 'saved')
    plt.close(fig)


if __name__ == '__main__':
    np.random.seed(1)  # 1 for nature, 2 for forecast ensembles
    n_ens = 4

    maindir = '/jetfs/home/lkugler/wrf_profiles/TEST/'
    dir_out = maindir + '/wrf/ens/2022-03-31/'
    os.makedirs(dir_out, exist_ok=True)

    save_csv = False
    plot = False
    prefix = 'raso.nat'

    for iens in range(1, n_ens+1):
        print('iens', iens)

        z, p, t, rh, ws, wd = load_reference_profile()

        ws, wd = specify_wind_profile(p)  # overwrites existing

        t, rh, u, v = perturb_profile(z, t, rh, ws, wd)
        
        r = calc_mixing_ratio(p, t, rh)  # use RH column as input moisture profile

        f_out = prefix+'.'+str(iens).zfill(3)+'.wrfprof'

        # use the mixing ratio column as input moisture profile
        # r = data[:, 5]/1000.
        kappa_moist = Rd/cpd*(1+r*Rv/Rd)/(1+r*cpv/cpd)
        pot_tmp = t*(1000/p)**kappa_moist  # potential temperature

        write_WRF_format(z, p[0], pot_tmp, r, u, v, f_out=dir_out+'/'+f_out)

        if save_csv:
            csvname = dir_out+'/'+'.'.join(f_out.split('.')[:-1])+'.csv'
            df = pd.DataFrame(data={'z': z, 'p': p, 'T': t, 'r': r*1000, 'U': u, 'V': v})
            df.to_csv(csvname)
            print(csvname, 'saved.')

        if plot:
            hodograph(z, u, v, 
                f_out=dir_out+'/'+'.'.join(f_out.split('.')[:-1])+'_hodograph.png')
            
            profile(p, t, r, u, v, 
                f_out=dir_out+'/'+'.'.join(f_out.split('.')[:-1])+'_skewt.png')
