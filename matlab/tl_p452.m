function Lb = tl_p452(f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp)
%tl_p452 basic transmission loss according to ITU-R P.452-18
%   Lb = tl_p452(f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp)
%
%   This is the MAIN function that computes the basic transmission loss not exceeded for p% of time
%   as defined in ITU-R P.452-18 (Section 4.5) for clear-air conditions. Other functions called from
%   this function are in ./private/ subfolder.
%
%     Input parameters:
%     f       -   Frequency (GHz)
%     p       -   Required time percentage for which the calculated basic
%                 transmission loss is not exceeded
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%                 above mean sea level. Both vectors contain n+1 profile points
%     g       -   vector of clutter + terrain profile heights gi along the path gi = hi + Ri (masl) 
%                 where Ri is the (representative) clutter height 
%     zone    -   vector of zone types: Coastal land (1), Inland (2) or Sea (3)
%     htg     -   Tx Antenna center heigth above ground level (m)
%     hrg     -   Rx Antenna center heigth above ground level (m)
%     phit_e  -   Tx Longitude (degrees)
%     phit_n  -   Tx Latitude  (degrees)
%     phir_e  -   Rx Longitude (degrees)
%     phir_n  -   Rx Latitude  (degrees)
%     Gt, Gr  -   Antenna gain in the direction of the horizon along the
%                 great-circle interference path (dBi)
%     pol     -   polarization of the signal (1) horizontal, (2) vertical
%     dct     -   Distance over land from the transmit and receive
%     dcr         antennas to the coast along the great-circle interference path (km).
%                 Set to zero for a terminal on a ship or sea platform
%     press   -   Dry air pressure (hPa)
%     temp    -   Air temperature (degrees C)
%
%     Output parameters:
%     Lb     -   basic  transmission loss according to ITU-R P.452-18
%
%     Example:
%     Lb = tl_p452(f, p, d, h, g, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp)


%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version
%     v1    10MAR16     Ivica Stevanovic, OFCOM         Introduced clutter losses
%     v2    15NOV16     Ivica Stevanovic, OFCOM         Corrected definition for Fj
%     v3    16NOV16     Ivica Stevanovic, OFCOM         Introduced the new version of ITU-R P.676-11 for gasseous attenuation calculation
%                                                       Corrected bug (epsr <-> sigma) in dl_se_ft
%                                                       Added a machine precision limit to avoid division by zero in tl_anomalous 
%                                                       Added a check for the number of points in the path profile 
%                                                       Added a check for the clutter loss nominal distances in cl_loss 
%     v4    10JAN17     Ivica Stevanovic, OFCOM         Corrected the Input parameters definition pol = 1 (horizontal), pol = 2 (vertical)
%     v5    13FEB17     Ivica Stevanovic, OFCOM         included lower limit for alpha and upper limit for mu2 in tl_anomalous
%     v6    05JUN20     Ivica Stevanovic, OFCOM         Introduced 3D distance in Free-space calculation
%                                                       Introduced a new computationally efficient version of find_intervals.m to align with ITU-R P.1812-5
%     v7    13JUL21     Ivica Stevanovic, OFCOM         Renamed subfolder "src" into "private" which is automatically in the MATLAB search path
%                                                       (as suggested by K. Konstantinou, Ofcom UK)   
%     v8    08OCT21     Ivica Stevanovic, OFCOM         Ensured that the variable "series" is a row vector in find_intervals.m
%     v9    24MAR22     Ivica Stevanovic, OFCOM         Introduced path center latitude as input argument (instead of Tx/Rx latitudes) 
%     v10   15MAR23     Ivica Stevanovic, OFCOM         Introduced troposcatter prediction model according to 3M/364 Annex 2
%     v11   25APR23     Ivica Stevanovic, OFCOM         Introduced more efficient interpolation and aligned with the rest of the Recommendation
%     v12   15NOV23     Ivica Stevanovic, OFCOM         Aligned with ITU-R P.452-18 (distributed clutter model) 

% MATLAB Version '9.12.0.1975300 (R2022a) Update 3' used in development of this code
%
% The Software is provided "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO, THE WARRANTIES OF MERCHANTABILITY, FITNESS 
% FOR A PARTICULAR PURPOSE AND NON-INFRINGMENT OF INTELLECTUAL PROPERTY RIGHTS 
% 
% Neither the Software Copyright Holder (or its affiliates) nor the ITU 
% shall be held liable in any event for any damages whatsoever
% (including, without limitation, damages for loss of profits, business 
% interruption, loss of information, or any other pecuniary loss)
% arising out of or related to the use of or inability to use the Software.
%
% THE AUTHOR(S) AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
%
% This function calls other functions that are placed in the ./private folder


% verify input argument values and limits
check_limit(f, 0.1, 50.0, 'f [GHz]');
check_limit(p, 0.001, 50, 'p [%]');
% check_limit(phi_path, -90, 90, 'phi_path [deg]'); todo: check longitudes and latitudes
check_limit(dct, 0, inf, 'dct [km]');
check_limit(dcr, 0, inf, 'dcr [km]');
check_value(pol, [1, 2], 'Polarization (pol) ');
check_value(zone, [1, 2, 3], 'Radio-climatic zone (zone) ');

if (d(1) > 0)
   error('d(0)  must be equal to zero.'); 
end

% Apply the condition in Step 4: Radio profile 
% gi is the terrain height in metres above sea level for all the points at a distance from transmitter or receiver less than 50 m.

kk = find(d < 50/1000);
if (~isempty(kk))
    g(kk) = h(kk);
end

endVal = d(end) - 50/1000;
kk = find(d > endVal);
if (~isempty(kk))
    g(kk) = h(kk);
end

% Compute the path profile parameters

% Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
zone_r = 12;
dtm = longest_cont_dist(d, zone, zone_r);

% Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
zone_r = 2;
dlm = longest_cont_dist(d, zone, zone_r);

% Calculate the longitude and latitude of the mid-point of the path, phim_e,
% and phim_n for dpnt = 0.5dt
Re = 6371;
dpnt = 0.5*(d(end)-d(1));
[phim_e, phim_n, bt2r, dgc] = great_circle_path(phir_e, phit_e, phir_n, phit_n, Re, dpnt);


% Find radio-refractivity lapse rate dN 
% using the digital maps at phim_e (lon), phim_n (lat) - as a bilinear interpolation

DN = get_interp2('DN50',phim_e,phim_n);
N0 = get_interp2('N050',phim_e,phim_n);


% Compute b0
b0 = beta0(phim_n, dtm, dlm);

[ae, ab] = earth_rad_eff(DN);

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);


[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Effective Earth curvature Ce (km^-1)

Ce = 1/ae;


% Find the intermediate profile point with the highest slope of the line
% from the transmitter to the point

if length(d)<4
    error('tl_p452: path profile requires at least 4 points.');
end

% Section 4.5: The overall prediction

di = d(2:end-1);
hi = h(2:end-1);

Stim = max((hi + 500*Ce*di.*(dtot - di) - hts)./di );           % Use hi instead of gi in Eq (14)

% Calculate the slope of the line from transmitter to receiver assuming a
% LoS path

Str = (hrs - hts)/dtot;                                         % Eq (15)

% Calculate an interpolation factor Fj to take account of the path angular
% distance (58)

THETA = 0.3;
KSI = 0.8;

Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (Stim-Str)/THETA) );

% Calculate an interpolation factor, Fk, to take account of the great
% circle path distance:

dsw = 20;
kappa = 0.5;

Fk = 1.0 - 0.5*( 1.0 + tanh(3.0 * kappa * (dtot-dsw)/dsw) ); % eq (59)

% modified with 3-D path for free-space computation, 
d3D = sqrt(dtot*dtot + ((hts-hrs)/1000.0).^2);

[Lbfsg, Lb0p, Lb0b] = pl_los(d3D, f, p, b0, omega, temp, press, dlt, dlr);


[ Ldp, Ld50 ] = dl_p( d, g, hts, hrs, hstd, hsrd, f, omega, p, b0, DN );

% The median basic transmission loss associated with diffraction Eq (43)

Lbd50 = Lbfsg + Ld50;

% The basic tranmission loss associated with diffraction not exceeded for
% p% time Eq (44)

Lbd = Lb0p + Ldp;

% A notional minimum basic transmission loss associated with LoS
% propagation and over-sea sub-path diffraction

Lminb0p = Lb0p + (1-omega)*Ldp;

if p >= b0
    
   Fi = inv_cum_norm(p/100)/inv_cum_norm(b0/100);   %eq (41a)
   
   Lminb0p = Lbd50 + (Lb0b + (1-omega)*Ldp - Lbd50)*Fi;   % eq (60)
   
end

% Calculate a notional minimum basic transmission loss associated with LoS
% and transhorizon signal enhancements

eta = 2.5;

Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, ae, b0);

Lminbap = eta*log(exp(Lba/eta) + exp(Lb0p/eta));    % eq (61)

% Calculate a notional basic transmission loss associated with diffraction
% and LoS or ducting/layer reflection enhancements

Lbda = Lbd;

if Lminbap <= Lbd
   Lbda = Lminbap + (Lbd-Lminbap)*Fk; 
end

% Calculate a modified basic transmission loss, which takes diffraction and
% LoS or ducting/layer-reflection enhancements into account

Lbam = Lbda + (Lminb0p - Lbda)*Fj;   % eq (63)

% Calculate the basic transmission loss due to troposcatter not exceeded for any time percantage p

Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr );

% Calculate the final transmission loss not exceeded for p% time

Lb_pol = -5*log10(10.^(-0.2*Lbs) + 10.^(-0.2*Lbam));  % eq (64)

Lb = Lb_pol(pol);


return
end