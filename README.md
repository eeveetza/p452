# MATLAB/Octave Implementation of Recommendation ITU-R P.452

[![DOI](https://zenodo.org/badge/459621140.svg)](https://zenodo.org/badge/latestdoi/459621140)


This code repository contains a MATLAB/Octave software implementation of Recommendation [ITU-R P.452-18](https://www.itu.int/rec/R-REC-P.452/en) with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 100 MHz. 

This implementation corresponds to the reference version approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).

<!--This development version implements the clutter loss model along the path profile.

This is development code that does not necessarily corresond to the reference version approved by ITU-R Working Party 3M and published by Study Group 3 on [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx).
-->

The following table describes the structure of the folder `./matlab/`.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_p452.m`                | MATLAB function implementing Recommendation ITU-R P.452-18         |
|`initiate_digital_maps.m`| MATLAB script that processes the ITU-R maps and generates the necessary functions. It needs to be run prior to using this software implementation. For details, see [Integrating ITU Digital Products](#integrating-itu-digital-products). |
|`validate_p452.m`          | MATLAB script used to validate the implementation of this Recommendation in `tl_p452.m` against the reference data in `validation_examples`.  It works in both MATLAB and Octave on Windows and MacOS.           |
|`./validation_examples/`    | Folder containing a non-exhaustive set of validation examples, in a form of .csv files, for different terrain profiles, clutter height profiles, frequencies, time-probabilities, etc. They include intermediate and final results of the calculations performed within P.452-18 with the aim of facilitating testing and validation, as well as comparison between different software implementations. |
|`./private/`   |  Folder containing all the MATLAB routines necessary for the implementation of the propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-11 (computing the specific attenuation due to dry air and water vapor by means of a summation of individual resonance lines from oxigen and water vapor). This folder contains test functions used to verify the current implementation of the model. It also contains several files with path profile data used in testing.|
|`P452.m`  `P452.fig`                | Graphical User Interface defined in those files can be opened by invoking the command `>> P452` in the MATLAB command window. Works only in MATLAB and not in Octave.       |
|`test_example.mat`                | An example of simulation data file that can be opened by the GUI. |


## Integrating ITU Digital Products

This software uses ITU digital products that are integral part of Recommendations. These products must not be reproduced or distributed without explicit written permission from the ITU.

### Setup Instructions

1. **Download and extract the required maps** to `./private/maps`:

   - From [ITU-R P.452-18](https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.452-18-202310-I!!ZIP-E.zip):
     - `N050.TXT`
     - `DN50.TXT`

2. **Run the script** `initiate_digital_maps.m` to generate the necessary functions for retrieving and interpolating data from from the maps.

### Notes

- Ensure all files are placed in `./private/maps` before running the script.
- The script processes the maps, which are critical for the software’s functionality.
- The resulting `*.m` files are placed in the folder `./private`.

## Function Call

The function `tl_p452` can be called by invoking the required input arguments:
~~~
Lb = tl_p452(f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp)
~~~


## Required input arguments of function `tl_p452`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `f`               | scalar double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `p         `      | scalar double | %     | 0.001 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | array double | km    |  0 < `max(d)` ≤ ~10000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | array double | m (asl)   |   | Terrain profile heights |
| `g`          | array double | m (asl)   |  | Clutter + Terrain profile heights   |
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



 
## Outputs ##

| Variable   | Type   | Units | Description |
|------------|--------|-------|-------------|
| `Lb`    | double | dB    | Basic transmission loss |


## Software Versions
The code was tested and runs on:
* MATLAB version 2022a 
* Octave version 6.1.0

## References

* [Recommendation ITU-R P.452](https://www.itu.int/rec/R-REC-P.452/en)

* [ITU-R SG 3 Software, Data, and Validation Web Page](https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx)
