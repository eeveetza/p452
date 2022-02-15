P452 Version 17.0_PDR (27.10.21)

GENERAL NOTES
----------------
- This MATLAB/Octave program computes the basic transmission loss according to PDR of Recommendation ITU-R P.452-16 that clutter along the path.

- The folder contains the following:


  1) tl_p452.m is the main function implementing PDR of Recommendation ITU-R P.452-16. This function calls
     other matlab routines defined in subfolder ./private/
     To know how to use the function, read the header of the function or type
     >> help tl_p452_pdr 

  2) Subfolder ./private/ with all the MATLAB routines necessary for the implementation of the
     propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-11
     (computing the specific attenuation due to dry air and water vapor by means of a summation
     of individual resonance lines from oxigen and water vapor).


All the scripts (except for the Graphical User Interface) work in Octave (versions 6.0 and above).

UPDATES AND FIXES
-----------------
Version 17.0_PDR (27.10.21)
    - Introduced clutter along the path as described in PDR P.452-16
    - Aligned pl_los.m with the ITU-R P.452-17 and free-space propagation model from ITU-R P.525-4
    - Updated validation examples to align with the changes in the factor 92.4 dB in free-space basic transmission loss
      Currently all the validation examples are without clutter profile - this needs tbd

Version 17.0 (08.10.21)
    - Renamed subfolder "src" into "private" which is automatically in the MATLAB search path
    - Modified the starting point in P452.m when computing transmission loss vs distance to make sure there are at least three points 
      between clutter at the Tx and Rx sides
    - Ensured that the variable series is a row vector in find_intervals.m

Version 16.3 (05.06.20)
    - Introduced 3D distance for free-space basic transmission loss (to be included in the next revision of Recommendation ITU-R P.452-16)
    - Introduced a new computationally efficient version of find_intervals.m to align with P.1812
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


