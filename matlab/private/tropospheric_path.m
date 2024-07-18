function [d_tcv, phi_cve, phi_cvn] = tropospheric_path(dt, hts, hrs, theta_e, theta_tpos, theta_rpos, phi_re, phi_te, phi_rn, phi_tn, Re)
%trophospheric path segments according to ITU-R P.2001-4
% [d_tcv, phi_cve, phi_cvn] = tropospheric_path(d, h, hts, theta_e, theta_tpos, theta_rpos, phi_re, phi_te, phi_rn, phi_tn, Re)
% This function computes tropospheric path segments as described in Section
% 3.9 of Recommendation ITU-R P.2001-4
%
% Input parameters:
% dt        -   Path length (km)
% hts, hrs  -   Tx/Rx antenna heights above means sea level (m)
% theta_e   -   Angle subtended by d km at the center of a sphere of effective earth radius (rad)
% theta_tpos-   Interfering antenna horizon elevation angle limited to be positive (mrad)
% theta_rpos-   Interfered-with antenna horizon elevation angle limited to be positive (mrad)
%               hts = htg + h(1)
% phi_re    -   Receiver longitude, positive to east (deg)
% phi_te    -   Transmitter longitude, positive to east (deg)
% phi_rn    -   Receiver latitude, positive to north (deg)
% phi_tn    -   Transmitter latitude, positive to north (deg)
% Re        -   Average Earth radius (km)
%
% Output parameters:
% d_tcv     -   Horizontal path length from transmitter to common volume (km)
% phi_cve   -   Longitude of the common volume 
% phi_cvn   -   Latitude of the common volume 

%
% Rev   Date        Author                          Description
% -------------------------------------------------------------------------------
% v0    13JUL16     Ivica Stevanovic, OFCOM         Initial version 
% v1    15MAR23     Ivica Stevanovic, OFCOM         Modified to tailor to P.452 PDR troposcatter

% Horizontal path length from transmitter to common volume (3.9.1a)

d_tcv = ( dt * tan(0.001*theta_rpos + 0.5* theta_e) - 0.001*(hts-hrs) ) / ...
        ( tan(0.001*theta_tpos + 0.5*theta_e) + tan(0.001*theta_rpos + 0.5* theta_e));

% Limit d_tcv such that 0 <= dtcv <= dt

if d_tcv < 0
    d_tcv = 0;
end
if d_tcv > dt
    d_tcv = dt;
end


% Calculate the longitude and latitude of the common volume from the
% transmitter and receiver longitudes and latitudes using the great circle
% path method of Attachment H by seting d_pnt = d_tcv

[phi_cve, phi_cvn, bt2r, dgc] = great_circle_path(phi_re, phi_te, phi_rn, phi_tn, Re, d_tcv);


return
end