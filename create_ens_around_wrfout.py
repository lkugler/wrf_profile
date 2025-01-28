import numpy as np
import xarray as xr
from create_wrf_profiles import write_WRF_format, perturb
from metpy import calc as mpcalc
from metpy.units import units
np.random.seed(1) 

def generate_ensemble(z, t, rh, u, v,
                      f_out='./test.txt'):

    t += perturb(0., .25, n_rand, z) * units.K
    rh += rh*perturb(0, .05, n_rand, z)

    u += perturb(0., .25, n_rand, z)
    v += perturb(0., .25, n_rand, z)

    rh[rh>1] = 1
    rh[rh<0] = 0

    # convert to mixing ratio
    r = mpcalc.mixing_ratio_from_relative_humidity(p, t, rh)
    
    write_WRF_format(z, psfc/100, t, r, u, v, f_out = f_out)

##################################################

# open the nature run file
f_wrfout = '/jetfs/home/lkugler/data/sim_archive/nat_250m_1600x1600x100/2008-07-30_08:00_/1/wrfout_d01_2008-07-30_10_00_00'

with xr.open_dataset(f_wrfout) as ds_in:
    # get surface pressure
    psfc = (ds_in['PSFC']).mean(dim=('south_north', 'west_east')).squeeze()

    # get PHB, PH, T, QVAPOR, U, V
    z = (ds_in['PHB'] + ds_in['PH']).mean(dim=('south_north', 'west_east')).squeeze()/9.81
    z0 = z.values[0]
    z = (z.values[1:] + z.values[:-1])/2
    p = (ds_in['P']+ ds_in['PB']).mean(dim=('south_north', 'west_east')).squeeze()

    t = (ds_in['T']).mean(dim=('south_north', 'west_east')).squeeze() + 300
    r = (ds_in['QVAPOR']).mean(dim=('south_north', 'west_east')).squeeze()
    u = (ds_in['U']).mean(dim=('south_north', 'west_east_stag')).squeeze()
    v = (ds_in['V']).mean(dim=('south_north_stag', 'west_east')).squeeze()

# write to text file
dir_out = '/jetfs/home/lkugler/data/sim_archive/nat_250m_1600x1600x100/2008-07-30_08:00_/1/'
f_out = dir_out+'/input_sounding_derived.txt'
write_WRF_format(z, psfc/100, t, r, u, v, f_out = f_out)

p = p*units.pascal
t = t*units.K
z = z*units.meter

## perturb the profile to create ensemble members
n_rand = 10  # number_of_random_numbers
n_ens = 40
    
# convert to RH
rh = mpcalc.relative_humidity_from_mixing_ratio(p, t, r)

f_out = '/jetfs/home/lkugler/wrf_profiles/data/2025-01-28/nat_250m_ens<iens>.wrfprof'
    
for iens in range(1, n_ens+1):
    print('iens', iens)
    # perturbations
    
    generate_ensemble(z, t, rh, u, v, 
                      f_out=f_out.replace('<iens>', str(iens).zfill(3)))

