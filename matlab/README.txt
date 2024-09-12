
P452 Version 18.1 (12.09.24)

GENERAL NOTES
----------------
- This MATLAB/Octave program computes the basic transmission loss according to Recommendation ITU-R P.452-18.

- The folder contains the following:


  1) tl_p452.m is the main function implementing Recommendation ITU-R P.452-18. This function calls
     other matlab routines defined in subfolder ./private/
     To know how to use the function, read the header of the function or type
     >> help tl_p452 

  2) initiate_digital_maps.m is the MATLAB script that processes the ITU-R maps and generates the necessary functions. 
     This software uses ITU digital products that are integral part of Recommendations. 
     These products must not be reproduced or distributed without explicit written permission from the ITU.

     a) Download and extract the required maps to `./private/maps` 
        From https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.452-18-202310-I!!ZIP-E.zip:
        - N050.TXT
        - DN50.TXT

     b) Run the script initiate_digital_maps.m to generate the necessary functions for retrieving and interpolating data. 
        The resulting `*.m` files are placed in the folder `./private`.
   
  3) Subfolder ./private/ with all the MATLAB routines necessary for the implementation of the
     propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-11
     (computing the specific attenuation due to dry air and water vapor by means of a summation
     of individual resonance lines from oxigen and water vapor).

  4) Graphical User Interface defined in files P452.m and P452.fig (only MATLAB)
     The interface can be opened by invoking the following command in the MATLAB command window 
     >> P452

  5) test_example.mat is an example of simulation data file that can be opened by the GUI
  
  6) Subfolder ./private/ also contains test functions used to verify the current implementation of the model.
     including several files with path profile data used in testing.

  7) Folder validation_examples containing the path profiles and reference data (as .csv files) to test the software implementation using validate_p452.m,
      which runs on both MATLAB and Octave, Windows and MacOS.

All the scripts (except for the Graphical User Interface) work in Octave (versions 6.0 and above).

UPDATES AND FIXES
-----------------
Version 18.1 (12.09.24)
    - Introduced functions and workflow for integrating ITU-R maps

Version 18.0 (16.05.24)
    - Aligned with Rec ITU-R P.452-18 (distributed clutter model)
    - Included separate input arguments for Tx/Rx longitude and latitude
    - Removed input arguments for surface refractivity and refractity gradient
    - Updated the MATLAB P452 GUI to account for the above changes
    - Updated the existing validation examples and added new validation examples

Version 17.1 (25.04.23)
    - Corrected an indexing issue in path_fraction.m and longest_cont_dist.m
    - Corrected various issues in pw2p.m

Version 17.0 (22.05.22)
    - Renamed subfolder "src" into "private" which is automatically in the MATLAB search path
    - Modified the starting point in P452.m when computing transmission loss vs distance to make sure there are at least three points 
      between clutter at the Tx and Rx sides
    - Ensured that the variable series is a row vector in find_intervals.m
    - Updated validation examples to align with the changes in the factor 92.4 dB in free-space basic transmission loss
    - Introduced validation examples in .csv instead of .xlsx 
    - Included path center latitude as input argument instead of Tx/Rx latitudes
    - Use 2.998e8 m/s as speed of light as per ITU-R P.2001-4 (instead of 3e8 m/s)
    - Updated P452 GUI to account for the input argument changes

Version 16.3 (05.06.20)
    - Introduced 3D distance for free-space basic transmission loss (to be included in the next revision of Recommendation ITU-R P.452-16)
    - Introduced a new computationally efficient version of find_intervals.m to align with P.1812-5
    - Introduced Octave specific code for reading/writing excel files in validation examples that now work in Octave as well

Version 2.15.02.17
    - Corrected the Input parameters definition pol = 1(horizontal), pol = 2 (vertical)
    - Included lower limit for alpha and upper limit for mu2 in tl_anomalous

Version 2.25.11.16
    - Introduced the new version of ITU-R P.676-11 for gaseous attenuation calculation
    - Corrected definition for Fj in tl_p452
    - Corrected bug (epsr <-> sigma) in dl_se_ft
    - Added a machine precision limit to avoid division by zero in tl_anomalous
    - Added a check for the number of points in the path profile
    - Added a check for the clutter loss nominal distances in cl_loss

Version 1.16.06.16
     - Initial implementation including clutter losses


License and copyright notice

Swiss Federal Office of Communications OFCOM (hereinafter the "Software Copyright Holder") makes the accompanying software 
(hereinafter the "Software") available free from copyright restriction. 

The Software Copyright Holder represents and warrants that to the best of its knowledge, 
it has the necessary copyright rights to waive all of the copyright rights as permissible under national law in the Software 
such that the Software can be used by implementers without further licensing concerns. 

No patent licence is granted, nor is a patent licensing commitment made, by implication, estoppel or otherwise. 

Disclaimer: Other than as expressly provided herein, 

(1) the Software is provided "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO, 
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGMENT OF INTELLECTUAL PROPERTY RIGHTS and 

(2) neither the Software Copyright Holder (or its affiliates) nor the ITU shall be held liable in any event for any damages whatsoever 
(including, without limitation, damages for loss of profits, business interruption, loss of information, or any other pecuniary loss) 
arising out of or related to the use of or inability to use the Software.


