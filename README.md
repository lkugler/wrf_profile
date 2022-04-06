# wrf_profile
Create input_sounding files for idealized WRF runs

1) Specify reference profile (e.g. by hand, https://github.com/lkugler/wrf_profile/blob/master/data/LMU/improved/raso.v2)
2) Run `create_wrf_profiles.py` (modify as you like)

#### Legacy code: alternative data source (Uni Wyoming) 
Use `convert_uwyo_WRF.py` with text file row-column data.
Remove trailing lines which contain missing values.
