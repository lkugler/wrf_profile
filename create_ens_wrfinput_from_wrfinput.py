import shutil
import numpy as np
import xarray as xr
import netCDF4 as nc
from metpy.units import units
from create_ens_profile_from_wrfinput import perturb

np.random.seed(1) 


def main(f_wrfinput_nature, f_wrfinput_ens, n_ens, number_of_random_numbers, f_archive_script):

    # extract the nature run profile (spatial mean)
    with xr.open_dataset(f_wrfinput_nature) as ds_in:
        # get surface pressure
        #psfc = (ds_in['PSFC']).mean(dim=('south_north', 'west_east')).squeeze()

        # get PHB, PH, T, QVAPOR, U, V
        z = (ds_in['PHB'] + ds_in['PH']).mean(dim=('south_north', 'west_east')).squeeze()/9.81 * units.meter
        #z0 = z.values[0]
        z = (z.values[1:] + z.values[:-1])/2
        #p = (ds_in['P']+ ds_in['PB']).mean(dim=('south_north', 'west_east')).squeeze() * units.pascal

        #theta = (ds_in['THM']).mean(dim=('south_north', 'west_east')).squeeze() *units.K
        r = (ds_in['QVAPOR']).mean(dim=('south_north', 'west_east')).squeeze() * units('kg/kg')
        #u = (ds_in['U']).mean(dim=('south_north', 'west_east_stag')).squeeze() * units('m/s')
        #v = (ds_in['V']).mean(dim=('south_north_stag', 'west_east')).squeeze() * units('m/s')


    nlevels = len(z)

    # define altitudes at which perturbation is maximum (at which to generate random numbers)
    z_top_perturb = (z.min() + 15000)
    z_coarse = np.linspace(z.min(), z_top_perturb, number_of_random_numbers)
    print('maximum of perturbation amplitude at this heights:', z_coarse)

    mask_perturb = z < z_top_perturb
    z_interp = z[mask_perturb]
    print('z of interpolated perturbations:', z_interp)

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
        theta_ens[iens, mask_perturb] = perturb(0., 1.0, z_coarse, z[mask_perturb]) * units.K # 1.0 K
        r_ens[iens, mask_perturb] = r[mask_perturb] * perturb(0, 0.01, z_coarse, z[mask_perturb]) # 0.05

        u_ens[iens, mask_perturb] = perturb(0., 0.3, z_coarse, z[mask_perturb]) * units('m/s') # 0.3
        v_ens[iens, mask_perturb] = perturb(0., 0.3, z_coarse, z[mask_perturb]) * units('m/s')  # 0.3

    # add the reversed perturbations
    theta_ens[n:, :] = -theta_ens[:n, :]
    r_ens[n:, :] = -r_ens[:n, :]
    u_ens[n:, :] = -u_ens[:n, :]
    v_ens[n:, :] = -v_ens[:n, :]

    # assert mean is zero along the ensemble dimension
    # assert that ensemble perturbations are unbiased
    assert np.allclose(theta_ens.mean(axis=0), 0)
    assert np.allclose(r_ens.mean(axis=0), 0)
    assert np.allclose(u_ens.mean(axis=0), 0)
    assert np.allclose(v_ens.mean(axis=0), 0)

    assert np.allclose(theta_ens[0,:], -theta_ens[n,:])
    assert np.allclose(theta_ens[1,:], -theta_ens[n+1,:])
    assert np.allclose(theta_ens[n-1,:], -theta_ens[-1,:])


    for iens in range(n_ens):
        print('iens', iens)
        
        # to get absolute values, add the profile perturbations to each of the ensemble member states 
        # (ensemble member states = PBL perturbed nature runs)
        f_in = f_wrfinput_ens.replace('<iens>', str(iens+1))
        
        # modifying this dataset
        with nc.Dataset(f_in, 'r+') as ds_ens:
        
            ds_ens.variables['THM'][:] += theta_ens[iens,:][np.newaxis, :, np.newaxis, np.newaxis]
            ds_ens.variables['QVAPOR'][:] += r_ens[iens,:][np.newaxis, :, np.newaxis, np.newaxis]
            ds_ens.variables['U'][:] += u_ens[iens,:][np.newaxis, :, np.newaxis, np.newaxis]
            ds_ens.variables['V'][:] += v_ens[iens,:][np.newaxis, :, np.newaxis, np.newaxis]
            
            print(f_in, 'modified')
            
    # for reproducibility, copy the script to the output directory
    shutil.copy(__file__, f_archive_script)
        

if __name__ == '__main__':
    """Create an ensemble of wrfinput files from a nature run wrfinput file.
    
    1) Run ideal.exe with the nature run input_profile, to get PBL perturbations of the nature run.
       If you want an ensemble, create one wrfinput_d01 file for each ensemble member using idea.exe.
       You can do this step by (future). 
    2) Run this script to add perturbations to each of the ensemble members' wrfinput files.
    
    Parameters:
    -----------
    f_wrfinput_nature : Path to nature run wrfinput file (created by ideal.exe)
    f_wrfinput_ens : Path to ensemble wrfinput file 
        (<iens> will be replaced by the ensemble number)
        Create these files in step 1. This script modifies the file.
    n_ens : Number of ensemble members to create
    
    Returns:
    --------
    None, f_wrfinput_ens will be overwritten. 
    """
    
    n_ens = 40
    number_of_random_numbers = 10

    # open the nature run file
    f_wrfinput_nature = '/jetfs/home/lkugler/data/sim_archive/nat_250m_1600x1600x100/2008-07-30_08:00/1/wrfinput_d01'
    
    # ensemble wrfinputs
    f_wrfinput_ens = '/jetfs/home/lkugler/data/run_WRF/exp_nat250m_noDA_2/<iens>/wrfinput_d01'
    
    f_archive_script = '/jetfs/home/lkugler/data/sim_archive/exp_nat250m_noDA_2/create_ens_perturbations.py'
    
    main(f_wrfinput_nature, f_wrfinput_ens, n_ens, number_of_random_numbers, f_archive_script)