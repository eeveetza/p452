# MATLAB/Octave Implementation of Preliminary Draft Recommendation (PDR) ITU-R P.452-17 with Distributed Clutter Loss Model and Corner Case Handling

This code repository contains a MATLAB/Octave software implementation of PDR [ITU-R P.452-17](https://www.itu.int/rec/R-REC-P.452/en) with a prediction procedure for the evaluation of interference between stations on the surface of the Earth at frequencies above about 0.1 GHz, as contained in Annex 8 of Working Party 3M Chairman's Report [3M/253 Annex 8](https://www.itu.int/dms_ties/itu-r/md/19/wp3m/c/R19-WP3M-C-0253!N08!MSW-E.docx), implementing clutter loss along the path on top of ITU-R P.452-17 (instead of terminal clutter loss). 
This code implements corner case handling.
This is a development version of the code and it is not approved by ITU-R Working Party 3M.


The following table describes the structure of the folder `./matlab/` containing the MATLAB/Octave implementation of Recommendation ITU-R P.452.

| File/Folder               | Description                                                         |
|----------------------------|---------------------------------------------------------------------|
|`tl_p452_pdr.m`                | MATLAB function implementing Recommendation ITU-R P.452-17_PDR          |
|`validate_p452_pdr.m`          | MATLAB script used to validate the implementation of this Recommendation in `tl_p452_pdr.m`.           |
|`./validation_examples/`    | Folder containing a non-exhaustive set of validation examples, in a form of .csv files, for different terrain profiles, frequencies, time-probabilities, etc. They include intermediate and final results of the calculations performed within P.452-17 with the aim of facilitating testing and validation, as well as comparison between different software implementations. Only terrain profile without clutter is used for testing at this moment.|
|`./private/`   |  Folder containing all the MATLAB routines necessary for the implementation of the propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-11 (computing the specific attenuation due to dry air and water vapor by means of a summation of individual resonance lines from oxigen and water vapor). This folder contains test functions used to verify the current implementation of the model. It also contains several files with path profile data used in testing.|




## Function Call

The function `tl_p452_pdr` can be called as follows:
~~~
Lb = tl_p452_pdr(f, p, d, h, zone, g, htg, hrg, phi_path, Gt, Gr, pol, dct, dcr, DN, N0, press, temp);
~~~


## Required input arguments of function `tl_p452`

| Variable          | Type   | Units | Limits       | Description  |
|-------------------|--------|-------|--------------|--------------|
| `f`               | scalar double | GHz   | ~0.1 ≤ `f` ≤ ~50 | Frequency   |
| `p         `      | scalar double | %     | 0.001 ≤ `p` ≤ 50 | Time percentage for which the calculated basic transmission loss is not exceeded |
| `d`               | array double | km    |  0 < `max(d)` ≤ ~10000 | Terrain profile distances (in the ascending order from the transmitter)|
| `h`          | array double | m (asl)   |   | Terrain profile heights |
| `zone`           | array int    |       | 1 - Coastal Land, 2 - Inland, 3 - Sea             |  Radio-climatic zone types |
| `g`           | array double    |       |              |  Clutter heights along the path |
| `htg`           | scalar double    | m      |           |  Tx antenna height above ground level |
| `hrg`           | scalar double    | m      |          |  Rx antenna height above ground level |
| `phi_path`           | scalar double    | deg      |   -90 ≤ `phi_path`  ≤ 90          |  Latitude of the path center between Tx and Rx stations |
| `Gt`  `Gr`           | scalar double  |   dBi    |           |  Tx/Rx antenna gain in the direction of the horizon towards along the great-circle interference path. |
| `pol`           | scalar int    |       |   `pol`  = 1, 2          |  Polarization of the signal: 1 - horizontal, 2 - vertical |
| `dct`           | scalar double    | km      |   `dct` ≥ 0          |  Distance over land from the Tx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `dcr`           | scalar double    | km      |   `dcr` ≥ 0          |  Distance over land from the Rx antenna to the coast along the great-circle interference path. To be set to zero for a terminal on a ship or sea platform.|
| `DN`            | scalar double    | N-units/km      | `DN`> 0           | The average radio-refractive index lapse-rate through the lowest 1 km of the atmosphere at the path-center. It can be derived from an appropriate map.  |
| `N0`           | scalar double    | N-units      |             | The sea-level surface refractivity at the path-centre. It can be derived from an appropriate map.|
| `press`           | scalar double    | hPa      |             | Dry air pressure.|
| `temp`           | scalar double    | deg C      |             | Air temperature.|



 
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
