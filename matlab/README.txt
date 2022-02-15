P452 Version 2 (15.02.17)

GENERAL NOTES
----------------
- This MATLAB program computes the basic transmission loss according to Recommendation ITU-R P.452-16.

- The folder contains the following:

  1) Graphical User Interface defined in files P452.m and P452.fig
     The interface can be opened by invoking the following command in the MATLAB command window
     >> P452

  2) tl_p452.m is the main function implementing Recommendation ITU-R P.452-16. This function calls
     other matlab routines defined in subfolder ./src/
     To know how to use the function, read the header of the function or type
     >> help tl_p452 

  3) Subfolder ./src/ with all the MATLAB routines necessary for the implementation of the
     propagation model, including the MATLAB implementation of Recommendation ITU-R P.676-10
     (computing the specific attenuation due to dry air and water vapor by means of a summation
     of individual resonance lines from oxigen and water vapor).

  4) Subfolder ./validation_examples with a non-exhaustive set of validation examples,
     in a form of Excel worksheets, for different terrain profiles, clutter heights, frequencies, time-probabilities, etc. 
     They include intermediate and final results of the calculations performed within P.452-16 
     with the aim of facilitating testing and validation, as well as comparison between different software implementations.
     

  5) validate_p452.m is a MATLAB script that validates the implementation of this recommendation as defined in tl_p452.m
     against the reference results obtained using Excel implementation of the same recommendation and reports success if
     the maximum deviation in final and intermediate results is less than tol = 1e-8.

  6) test_example.mat is an example of simulation data file that can be opened by the GUI
  
  7) Subfolder ./test/ with test functions used to verify the current implementation of the model.
     This subfolder also contains several files with path profile data used in testing.



UPDATES AND FIXES
-----------------
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



