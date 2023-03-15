# MATLAB/Octave Implementation of Recommendation ITU-R P.452


This code repository contains a MATLAB/Octave software implementation of Recommendation [ITU-R P.452-17](https://www.itu.int/rec/R-REC-P.452/en) with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz.  

This development version implements PDR from WP 3M Chairman's Report [3M/364 Annex 2](https://www.itu.int/dms_ties/itu-r/md/19/wp3m/c/R19-WP3M-C-0364!N02!MSW-E.docx) with troposcatter transmission loss prediction model from [ITU-R P.617](https://www.itu.int/rec/R-REC-P.617/en), which is being evaluated by CG 3K-3M-18.

This is development code that does not necessarily corresond to the reference version approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).

The following table describes the structure of the folder `./matlab/`.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_p452.m`                | MATLAB function implementing Recommendation ITU-R P.452-17          |
|`./C3_1_profiles/`   |             Folder containing terrain profiles and measurement files for Terrestrial trans-horizon links in DBSG3 (`CG-3M-2/DBSG3 Repository/Part II Terrestrial trans-horizon.../ III-01)`|
|`read_C3_1_profile.m`   |             MATLAB script for reading the terrain profile|
|`read_data_table_C3_1.m`   |             MATLAB script for reading the measurement data|
|`compute_btl_table_C3_1.m`   |             MATLAB script for computing the basic transmission loss according to ITU-R P.452-17 and the PDR on troposcatter for those paths from table C3_1 which have complete set of input parameters and computes the prediction errors of the two approaches in an Excel file `Results_Table_C3_1_P452.xls`|
|`validate_p452.m`          | MATLAB script used to validate the implementation of this Recommendation in `tl_p452.m` against the reference data in `validation_examples`.  It works in both MATLAB and Octave on Windows and MacOS.           |
|`./validation_examples/`    | Folder containing a non-exhaustive set of validation examples, in a form of .csv files, for different terrain profiles, clutter heights, frequencies, time-probabilities, etc. They include intermediate and final results of the calculations performed within P.452-17 with the aim of facilitating testing and validation, as well as comparison between different software implementations. |
|`./private/`   |  Folder containing all the MATLAB routines necessary for the implementation of the propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-11 (computing the specific attenuation due to dry air and water vapor by means of a summation of individual resonance lines from oxigen and water vapor). This folder contains test functions used to verify the current implementation of the model. It also contains several files with path profile data used in testing.|
|`P452.m`  `P452.fig`                | Graphical User Interface defined in those files can be opened by invoking the command `>> P452` in the MATLAB command window. Works only in MATLAB and not in Octave.       |
|`test_example.mat`                | An example of simulation data file that can be opened by the GUI. Works only in MATLAB and not in Octave.|




## Function Call

The function `tl_p452` can be called

1. by invoking only the required input arguments:
~~~
Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, pdr)
~~~

2. or by explicitly invoking all the input arguments (both required and optional):
~~~
Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, pdr, ...
             ha_t, ha_r, dk_t, dk_r);
~~~

## Required input arguments of function `tl_p452`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `f`               | scalar double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `p         `      | scalar double | %     | 0.001 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | array double | km    |  0 < `max(d)` ≤ ~10000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | array double | m (asl)   |   | Terrain profile heights |
| `zone`           | array int    |       | 1 - Coastal Land, 2 - Inland, 3 - Sea             |  Radio-climatic zone types |
| `htg`           | scalar double    | m      |           |  Tx antenna height above ground level |
| `hrg`           | scalar double    | m      |          |  Rx antenna height above ground level |
| `phit_e`           | scalar double    | deg      |     0 ≤ `phit_e`  ≤ 360          |  Tx longitude |
| `phit_n`           | scalar double    | deg      |     -90 ≤ `phit_n`  ≤ 90          |  Tx latitude |
| `phir_e`           | scalar double    | deg      |     0 ≤ `phir_e`  ≤ 360          |  Rx longitude |
| `phir_n`           | scalar double    | deg      |     -90 ≤ `phir_n`  ≤ 90          |  Rx latitude |
| `Gt`  `Gr`           | scalar double  |   dBi    |           |  Tx/Rx antenna gain in the direction of the horizon towards along the great-circle interference path. |
| `pol`           | scalar int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |
| `dct`           | scalar double    | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `dcr`           | scalar double    | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `press`           | scalar double    | hPa      |             | Dry air pressure.|
| `temp`           | scalar double    | deg C      |             | Air temperature.|
| `pdr`                 | scalar int    |        |   1, 0         |  0 - uses ITU-R P.452-17, 1 - uses PDR ITU-R P.452-17 troposcatter |

## Optional input arguments of function `tl_p452`
| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `ha_t`           | scalar double    | m      |             | Clutter nominal height at the Tx side |
| `ha_r`           | scalar double    | m      |             | Clutter nominal height at the Rx side |
| `dk_t`           | scalar double    | km      |             | Clutter nominal distance at the Tx side |
| `dk_r`           | scalar double    | km      |             | Clutter nominal distance at the Rx side |



 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | double | dB    | Basic transmission loss |


## Software Versions
The code was tested and runs on:
* MATLAB version 2022a (MS Windows)

## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)
