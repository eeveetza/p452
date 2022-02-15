P452 Version 16.3 (05.06.20)

GENERAL NOTES
----------------
- This MATLAB/Octave program computes the basic transmission loss according to Recommendation ITU-R P.452-16.

- The folder contains the following:


  1) tl_p452.m is the main function implementing Recommendation ITU-R P.452-16. This function calls
     other matlab routines defined in subfolder ./src/
     To know how to use the function, read the header of the function or type
     >> help tl_p452 

  2) Subfolder ./src/ with all the MATLAB routines necessary for the implementation of the
     propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-10
     (computing the specific attenuation due to dry air and water vapor by means of a summation
     of individual resonance lines from oxigen and water vapor).

  3) Subfolder ./validation_examples with a non-exhaustive set of validation examples,
     in a form of Excel worksheets, for different terrain profiles, clutter heights, frequencies, time-probabilities, etc. 
     They include intermediate and final results of the calculations performed within P.452-16 
     with the aim of facilitating testing and validation, as well as comparison between different software implementations.
     They are obtained using Excel implementation of Recommendation ITU-R P.452-16 as defined in ITUR_452_16.xlsm

  4) validate_p452.m is a MATLAB script that validates the implementation of this recommendation as defined in tl_p452.m
     against the reference results obtained using Excel implementation of the same recommendation and reports success if
     the maximum deviation in final and intermediate results is less than tol = 1e-4 dB.

  5) Graphical User Interface defined in files P452.m and P452.fig
     The interface can be opened by invoking the following command in the MATLAB command window 
     >> P452

  6) test_example.mat is an example of simulation data file that can be opened by the GUI
  
  7) Subfolder ./test/ with test functions used to verify the current implementation of the model.
     This subfolder also contains several files with path profile data used in testing.

All the scripts (except for the Graphical User Interface) should work in Octave as well.

UPDATES AND FIXES
-----------------
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
(hereinafter the “Software”) available free from copyright restriction. 

The Software Copyright Holder represents and warrants that to the best of its knowledge, 
it has the necessary copyright rights to waive all of the copyright rights as permissible under national law in the Software 
such that the Software can be used by implementers without further licensing concerns. 

No patent licence is granted, nor is a patent licensing commitment made, by implication, estoppel or otherwise. 

Disclaimer: Other than as expressly provided herein, 

(1) the Software is provided “AS IS” WITH NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO, 
THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGMENT OF INTELLECTUAL PROPERTY RIGHTS and 

(2) neither the Software Copyright Holder (or its affiliates) nor the ITU shall be held liable in any event for any damages whatsoever 
(including, without limitation, damages for loss of profits, business interruption, loss of information, or any other pecuniary loss) 
arising out of or related to the use of or inability to use the Software.


