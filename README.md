# wrf_profile
Create input_sounding files for idealized WRF runs

1) Specify reference profile (e.g. by hand, https://github.com/lkugler/wrf_profile/blob/master/data/LMU/improved/raso.v2)
2) Create perturbed initial condition profiles from reference profile with `create_perturbed_members.py` 
3) Transform to WRF input_profile format with `convert_LMU_modified_WRF.py` 

### alternative data source (Uni Wyoming)

```
f_in = 'wyoming_sounding.txt'  # excerpt of the html
f_out = 'input_sounding'  # to be used by ideal.exe of WRF


p, z, T, Td, r, dd, ff, pot_tmp = read_uwyo(f_in)

# convert to WRF variables
u = -ff * np.cos(dd/180*np.pi-np.pi/2)
v = -ff * np.sin(dd/180*np.pi+np.pi/2)

p_2m = p[0]  # surface pressure
write_wrfprof(p_2m, z, pot_tmp, r, u, v, f_out)
```

where `wyoming_sounding.txt` is like this:
(get this data from here: [weather.uwyo.edu](http://weather.uwyo.edu/upperair/europe.html))
or replace date and station here: [http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2020&MONTH=06&FROM=1912&TO=1912&STNM=10868](http://weather.uwyo.edu/cgi-bin/sounding?region=europe&TYPE=TEXT%3ALIST&YEAR=2020&MONTH=06&FROM=1912&TO=1912&STNM=10868)

```
10868 Muenchen-Oberschlssheim Observations at 12Z 19 Jun 2020

-----------------------------------------------------------------------------
   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
    hPa     m      C      C      %    g/kg    deg   knot     K      K      K 
-----------------------------------------------------------------------------
 1000.0    157                                                               
  961.0    492   12.8   12.2     96   9.37    215      2  289.2  315.7  290.9
  951.0    580   12.1   11.6     97   9.08    260     10  289.4  315.1  290.9
  925.0    812   10.2    9.9     98   8.34    270     21  289.7  313.5  291.2
  924.0    821   10.2    9.8     98   8.32    270     21  289.8  313.5  291.2
  855.0   1465    7.2    6.5     95   7.13    315     12  293.2  313.9  294.5
  850.0   1514    7.0    6.2     95   7.04    315     12  293.5  314.0  294.7
  845.0   1563    7.0    5.7     91   6.84    313     12  294.0  313.9  295.2
  826.0   1747    5.6    4.5     93   6.44    305     12  294.4  313.3  295.5
  760.0   2424    0.4    0.2     99   5.13    357      7  295.9  311.2  296.8
  710.0   2968   -1.8   -1.9    100   4.72     40      4  299.2  313.6  300.1
```
(stop at the last line without missing values)
