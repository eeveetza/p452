# MATLAB/Octave Implementation of Recommendation ITU-R P.452

[![DOI](https://zenodo.org/badge/459621140.svg)](https://zenodo.org/badge/latestdoi/459621140)


This code repository contains a MATLAB/Octave software implementation of Recommendation [ITU-R P.452-17](https://www.itu.int/rec/R-REC-P.452/en) with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz.  

This is code corresponds to the reference version approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).

The following table describes the structure of the folder `./matlab/` containing the MATLAB/Octave implementation of Recommendation ITU-R P.452.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_p452.m`                | MATLAB function implementing Recommendation ITU-R P.452-17          |
|`validate_p452.m`          | MATLAB script used to validate the implementation of this Recommendation in `tl_p452.m` against the reference data in `validation_examples`.  It works in both MATLAB and Octave on Windows and MacOS.           |
|`./validation_examples/`    | Folder containing a non-exhaustive set of validation examples, in a form of .csv files, for different terrain profiles, clutter heights, frequencies, time-probabilities, etc. They include intermediate and final results of the calculations performed within P.452-17 with the aim of facilitating testing and validation, as well as comparison between different software implementations. |
|`./private/`   |  Folder containing all the MATLAB routines necessary for the implementation of the propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-11 (computing the specific attenuation due to dry air and water vapor by means of a summation of individual resonance lines from oxigen and water vapor). This folder contains test functions used to verify the current implementation of the model. It also contains several files with path profile data used in testing.|
|`P452.m`  `P452.fig`                | Graphical User Interface defined in those files can be opened by invoking the command `>> P452` in the MATLAB command window. Works only in MATLAB and not in Octave.       |
|`test_example.mat`                | An example of simulation data file that can be opened by the GUI. Works only in MATLAB and not in Octave.|



## Function Call

The function `tl_p452` can be called

1. by invoking only the required input arguments:
~~~
Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_path, Gt, Gr, pol, dct, dcr, DN, N0, press, temp);
~~~

2. or by explicitly invoking all the input arguments (both required and optional):
~~~
Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_path, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, ...
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
| `phi_path`           | scalar double    | deg      |   -90 ≤ `phi_path`  ≤ 90          |  Latitude of path centre between Tx and Rx station |
| `Gt`  `Gr`           | scalar double  |   dBi    |           |  Tx/Rx antenna gain in the direction of the horizon towards along the great-circle interference path. |
| `pol`           | scalar int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |
| `dct`           | scalar double    | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `dcr`           | scalar double    | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `DN`            | scalar double    | N-units/km      | `DN`> 0           | The average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere at the path-center. It can be derived from an appropriate map.  |
| `N0`           | scalar double    | N-units      |             | The sea-level surface refractivity at the path-centre. It can be derived from an appropriate map.|
| `press`           | scalar double    | hPa      |             | Dry air pressure.|
| `temp`           | scalar double    | deg C      |             | Air temperature.|

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
* MATLAB versions 2017a and 2020a
* Octave version 6.1.0

## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)
