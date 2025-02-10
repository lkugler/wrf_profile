import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

from metpy import calc as mpcalc
from metpy.units import units

import matplotlib.pyplot as plt

np.random.seed(1) 
    
GRAVITY = 9.81 * units.meter / units.second**2
Rd = 287.05 * units.joule / units.kilogram / units.kelvin
Cp = 1005.7 * units.joule / units.kilogram / units.kelvin


def perturb(mu, sigma, z_coarse, z_interp):
    """get correlated perturbations on original grid
    
    Args:
        mu: float
        sigma: float
        number_of_random_numbers: int
        z: np.array
    
    Returns:
        np.array
    """
    number_of_random_numbers = len(z_coarse)

    # perturbations on coarse grid (uncorrelated)
    rand_coarse = np.random.normal(mu, sigma, number_of_random_numbers)
    for zi, ri in zip(z_coarse, rand_coarse):
        print('random pertubartion (height, value):', zi, ri)

    # interpolate in space to get autocorrelated perturbations
    f = interp1d(z_coarse, rand_coarse, kind='cubic')
    interpolated_randn  = f(z_interp)
    return interpolated_randn


def read_profile_from_wrfout(f_wrfout):
    with xr.open_dataset(f_wrfout) as ds_in:
        # get surface pressure
        psfc = (ds_in['PSFC']).mean(dim=('south_north', 'west_east')).squeeze()

        # get PHB, PH, T, QVAPOR, U, V
        z = (ds_in['PHB'] + ds_in['PH']).mean(dim=('south_north', 'west_east')).squeeze()/9.81
        z0 = z.values[0]
        z = (z.values[1:] + z.values[:-1])/2
        p = (ds_in['P']+ ds_in['PB']).mean(dim=('south_north', 'west_east')).squeeze()

        theta = (ds_in['T']).mean(dim=('south_north', 'west_east')).squeeze() + 300
        r = (ds_in['QVAPOR']).mean(dim=('south_north', 'west_east')).squeeze()
        u = (ds_in['U']).mean(dim=('south_north', 'west_east_stag')).squeeze()
        v = (ds_in['V']).mean(dim=('south_north_stag', 'west_east')).squeeze()
        
        
    p = p*units.pascal
    z = z*units.meter
    theta = theta*units.K
    return z, p, theta, r, u, v, psfc

def read_profile_from_input_sounding(f_in):
    """Read this text file into numpy arrays.
    
    Format is like this:
   962.00    301.22      13.79 
   531.83    301.22      13.79      -2.48       0.19
   617.76    301.07      13.64      -2.96       0.57
   704.26    301.02      13.79      -3.31       1.18
   791.48    301.18      14.73      -3.61       1.90
   879.44    301.39      14.63      -3.73       2.68
   
   """

    # read data with numpy
    column0 = np.loadtxt(f_in, usecols=(0))
    psfc = column0[0] *100 * units.pascal
    z = column0[1:] * units.meter

    column12 = np.loadtxt(f_in, usecols=(1, 2))
    theta = column12[:,0] * units.K
    r = column12[:,1] * units('g/kg')
    
    column34 = np.loadtxt(f_in, usecols=(3, 4), skiprows=1)
    u = column34[:,0] * units('m/s')
    v = column34[:,1] * units('m/s')
    
    # assert psfc is scalar
    assert psfc.size == 1
    
    # len(z) = len(theta)-1
    assert len(z) == len(theta)-1
    assert len(theta) == len(r)
    assert len(z) == len(u)
    assert len(z) == len(v)
    
    return psfc, z, theta, r, u, v
    

    
##################################################
if __name__ == '__main__':

    # open the nature run file
    # f_wrfout = '/jetfs/home/lkugler/data/sim_archive/nat_250m_1600x1600x100/2008-07-30_08:00_/1/wrfout_d01_2008-07-30_08_05_00'
    # z, p, theta, r, u, v, psfc = read_profile_from_wrfout(f_wrfout)

    # read state from input_sounding
    f_in = '/jetfs/home/lkugler/data/sim_archive/nat_250m_1600x1600x100/2008-07-30_08:00_/1/input_sounding'
    psfc, z, theta, r, u, v = read_profile_from_input_sounding(f_in)
    # first level of theta and r is duplicate (surface value)
    theta = theta[1:]
    r = r[1:]
    
    # define altitudes at which perturbation is maximum (at which to generate random numbers)
    number_of_random_numbers = 10
    z_top_perturb = z.min() + 15000 * units.meter
    z_coarse = np.linspace(z.min(), z_top_perturb, number_of_random_numbers)
    print('maximum of perturbation amplitude at this heights:', z_coarse)
    
    mask_perturb = z < z_top_perturb
    z_interp = z[mask_perturb]
    print('z of interpolated perturbations:', z_interp)

    # write to text file
    # dir_out = '/jetfs/home/lkugler/data/sim_archive/nat_250m_1600x1600x100/2008-07-30_08:00_/1/'
    # f_out = dir_out+'/input_sounding_derived.txt'
    # write_WRF_format(z, psfc/100, t, r, u, v, f_out = f_out)

    ## perturb the profile to create ensemble members
    n_rand = 10  # number_of_random_numbers
    n_ens = 40

    f_out = '/jetfs/home/lkugler/wrf_profiles/data/2025-02-03/nat_250m_ens<iens>.wrfprof'
        

    nlevels = len(z)

    # extend theta, r, u, v to 2d arrays with n_ens columns
    theta_ens = np.zeros((n_ens, nlevels))
    r_ens = np.zeros((n_ens, nlevels))
    u_ens = np.zeros((n_ens, nlevels))
    v_ens = np.zeros((n_ens, nlevels))
    
    # create n_ens/2 perturbations, assert n is integer
    assert n_ens % 2 == 0
    n = n_ens // 2

    for iens in range(n):
        print('iens', iens)

        # level 0 is the surface (should be equal to the first level)
        theta_ens[iens, mask_perturb] = 1 # perturb(0., 0., z_coarse, z[mask_perturb]) * units.K # 0.25
        r_ens[iens, mask_perturb] = r[mask_perturb] * perturb(0, 0., z_coarse, z[mask_perturb]) # 0.02

        u_ens[iens, mask_perturb] = perturb(0., 0., z_coarse, z[mask_perturb]) * units('m/s') # 0.25
        v_ens[iens, mask_perturb] = perturb(0., 0., z_coarse, z[mask_perturb]) * units('m/s')  # 0.25

    # add the reversed perturbations
    theta_ens[n:, :] = -theta_ens[:n, :]
    r_ens[n:, :] = -r_ens[:n, :]
    u_ens[n:, :] = -u_ens[:n, :]
    v_ens[n:, :] = -v_ens[:n, :]

    # assert mean is zero along the ensemble dimension
    assert np.allclose(theta_ens.mean(axis=0), 0)
    assert np.allclose(r_ens.mean(axis=0), 0)
    assert np.allclose(u_ens.mean(axis=0), 0)
    assert np.allclose(v_ens.mean(axis=0), 0)

    assert np.allclose(theta_ens[0,:], -theta_ens[n,:])
    assert np.allclose(theta_ens[1,:], -theta_ens[n+1,:])
    assert np.allclose(theta_ens[n-1,:], -theta_ens[-1,:])
    # rh[rh>.9] = .9
    # rh[rh<0] = 0

    #r = mpcalc.mixing_ratio_from_relative_humidity(p, t, rh)
    #theta = mpcalc.potential_temperature(p, t)
    
    
    # plot profiles of perturbations
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    for iens in range(n_ens):
        ax[0].plot(theta_ens[iens], z, color='r', alpha=0.5)
        ax[1].plot(r_ens[iens], z, color='g', alpha=0.5)
    
    ax[0].set_xlabel('∆theta [K]')
    ax[1].set_xlabel('∆r [g/kg]')
    fig.savefig('delta_profiles.png')
    
    
    # add reference states
    theta_ens = theta_ens*units.K + theta
    r_ens = r_ens*units('g/kg') + r
    u_ens = u_ens*units('m/s') + u
    v_ens = v_ens*units('m/s') + v
    
    
    # plot profiles of perturbations
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    for iens in range(n_ens):
        p = np.zeros(nlevels)
        t = np.zeros(nlevels) 
        # integrate to get pressure at each height
        p_below = psfc
        p_ref = 1e5*units.pascal
        t_lowest = theta_ens[iens, 0] * (p_below / p_ref)**(Rd / Cp)
        
        for i in range(1, nlevels-1):
            t[i] = (theta_ens[iens, i] * (p_below / p_ref)**(Rd / Cp)).magnitude * units.K
            p[i] = (p_below * np.exp(GRAVITY * (z[i-1] - z[i]) / Rd / t[i])).magnitude
            print('p[i]', p[i])
            p_below = p[i]
        
        t = mpcalc.temperature_from_potential_temperature(p, theta_ens[iens, :])
        ax[0].plot(t[iens, mask_perturb], z[mask_perturb], color='r', alpha=0.5)
        ax[1].plot(r_ens[iens, mask_perturb], z[mask_perturb], color='g', alpha=0.5)
    
    ax[0].set_xlabel('theta [K]')
    ax[1].set_xlabel('r [g/kg]')
    fig.savefig('profiles.png')
    
    
    from create_wrf_profiles import write_WRF_format
    
    for iens in range(n_ens):
        print('iens', iens)
        
        # add surface theta and r
        theta = theta_ens[iens, :]
        r = r_ens[iens, :]
        theta = np.concatenate([theta[:1], theta])
        r = np.concatenate([r[:1], r])
        
        u = u_ens[iens, :]
        v = v_ens[iens, :]
        
        write_WRF_format(z, psfc/100, theta, r, u, v, 
                         f_out = f_out.replace('<iens>', str(iens+1).zfill(3)))
