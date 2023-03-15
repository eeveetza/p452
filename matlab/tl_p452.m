function Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, pdr, varargin)
%tl_p452 basic transmission loss according to ITU-R P.452-17
%   Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, pdr, ha_t, ha_r, dk_t, dk_r )
%
%   This is the MAIN function that computes the basic transmission loss not exceeded for p% of time
%   as defined in ITU-R P.452-17 (Section 4.6). Other functions called from
%   this function are in ./private/ subfolder.
%
%     Input parameters:
%     f       -   Frequency (GHz)
%     p       -   Required time percentage for which the calculated basic
%                 transmission loss is not exceeded
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%                 above mean sea level. Both vectors contain n+1 profile points
%     zone    -   Zone type: Coastal land (1), Inland (2) or Sea (3)
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
%     pdr     -   pdr flag = 1, use PDR Troposcatter model
%     ha_t    -   Clutter nominal height (m) at the Tx side
%     ha_r    -   Clutter nominal height (m) at the Rx side
%     dk_t    -   Clutter nominal distance (km) at the Tx side
%     dk_r    -   Clutter nominal distance (km) at the Rx side
%
%     Output parameters:
%     Lb     -   basic  transmission loss according to ITU-R P.452-17
%
%     Example:
%     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, pdr)
%     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phit_e, phit_n, phir_e, phir_n, Gt, Gr, pol, dct, dcr, press, temp, pdr, ha_t, ha_r, dk_t, dk_r)

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
%     v4    10JAN17     Ivica Stevanovic, OFCOM         Corrected the Input parameters definition pol = 1(horizontal), pol = 2 (vertical) t
%     v5    13FEB17     Ivica Stevanovic, OFCOM         included lower limit for alpha and upper limit for mu2 in tl_anomalous
%     v6    05JUN20     Ivica Stevanovic, OFCOM         Introduced 3D distance in Free-space calculation
%                                                       Introduced a new computationally efficient version of find_intervals.m to align with ITU-R P.1812-5
%     v7    13JUL21     Ivica Stevanovic, OFCOM         Renamed subfolder "src" into "private" which is automatically in the MATLAB search path
%                                                       (as suggested by K. Konstantinou, Ofcom UK)   
%     v8    08OCT21     Ivica Stevanovic, OFCOM         Ensured that the variable "series" is a row vector in find_intervals.m
%     v9    24MAR22     Ivica Stevanovic, OFCOM         Introduced path center latitude as input argument (instead of Tx/Rx latitudes) 
%     v10   15MAR23     Ivica Stevanovic, OFCOM         Introduced troposcatter prediction model according to 3M/364 Annex 2
% MATLAB Version 9.7.0.1190202 (R2019b) used in development of this code
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


% Read the input arguments 

if nargin > 23    warning(strcat('tl_p452: Too many input arguments; The function requires at most 21',...
        'input arguments. Additional values ignored. Input values may be wrongly assigned.'));
end

if nargin <19 
    error('tl_p452: function requires at least 19 input parameters.');
end

ha_t = [];
ha_r = [];
dk_t = [];
dk_r = [];

narg = 20;


if nargin >=narg
    ha_t=varargin{1};
    narg = narg + 1;
    if nargin >=narg
        ha_r=varargin{2};
        narg = narg + 1;
        if nargin >=narg
            dk_t=varargin{3};
            narg = narg + 1;
            if nargin >=narg
                dk_r=varargin{4};
            end
        end
    end
end

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

DN50 = DigitalMaps_DN50();
N050 = DigitalMaps_N050();

latcnt = 90:-1.5:-90;               %Table 2.4.1
loncnt = 0:1.5:360;                 %Table 2.4.1  
[LON,LAT] = meshgrid(loncnt, latcnt);

% Map phicve (-180, 180) to loncnt (0,360);
phim_e1 = phim_e;
if phim_e1 < 0
    phim_e1 = phim_e + 360;
end

% Find radio-refractivity lapse rate dN 
% using the digital maps at phim_e (lon), phim_n (lat) - as a bilinear interpolation

DN = interp2(LON,LAT,DN50,phim_e1,phim_n);
N0 = interp2(LON,LAT,N050,phim_e1,phim_n);
clear DN50 N050

% Compute b0
b0 = beta0(phim_n, dtm, dlm);

[ae, ab] = earth_rad_eff(DN);

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

% Modify the path according to Section 4.5.4, Step 1 and compute clutter losses
% only if not isempty ha_t and ha_r

[dc, hc, zonec, htgc, hrgc, Aht, Ahr] = closs_corr(f, d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r);

d = dc;
h = hc;
zone = zonec;
htg = htgc;
hrg = hrgc;

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htgc;
hrs = h(end) + hrgc;

% Effective Earth curvature Ce (km^-1)

Ce = 1/ae;


% Find the intermediate profile point with the highest slope of the line
% from the transmitter to the point

if length(d)<4
    error('tl_p452: path profile requires at least 4 points.');
end

di = d(2:end-1);
hi = h(2:end-1);

Stim = max((hi + 500*Ce*di.*(dtot - di) - hts)./di );           % Eq (14)

% Calculate the slope of the line from transmitter to receiver assuming a
% LoS path

Str = (hrs - hts)/dtot;                                         % Eq (15)

% Calculate an interpolation factor Fj to take account of the path angular
% distance (58)

THETA = 0.3;
KSI = 0.8;

% changed the definition for Fj on 15DEC16.
%Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (theta-THETA)/THETA) )
Fj = 1.0 - 0.5*( 1.0 + tanh(3.0 * KSI * (Stim-Str)/THETA) );

% Calculate an interpolation factor, Fk, to take account of the great
% circle path distance:

dsw = 20;
kappa = 0.5;

Fk = 1.0 - 0.5*( 1.0 + tanh(3.0 * kappa * (dtot-dsw)/dsw) ); % eq (59)

% modified with 3-D path for free-space computation, 
d3D = sqrt(dtot*dtot + ((hts-hrs)/1000.0).^2);

%[Lbfsg, Lb0p, Lb0b] = pl_los(dtot, f, p, b0, omega, temp, press, dlt, dlr);
[Lbfsg, Lb0p, Lb0b] = pl_los(d3D, f, p, b0, omega, temp, press, dlt, dlr);


[ Ldp, Ld50 ] = dl_p( d, h, hts, hrs, hstd, hsrd, f, omega, p, b0, DN );

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

if (~pdr)

    Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr );

else

    % The path length expressed as the angle subtended by d km at the center of
    % a sphere of effective Earth radius ITU-R P.2001-4 (3.5.4)

    theta_e = dtot/ae; % radians

    % Calculate the horizon elevation angles limited such that they are positive

    theta_tpos = max(theta_t, 0);                   % Eq (3.7.11a) ITU-R P.2001-4
    theta_rpos = max(theta_r, 0);                   % Eq (3.7.11b) ITU-R P.2001-4

    [dt_cv, phi_cve, phi_cvn] = tropospheric_path(dtot, hts, hrs, theta_e, theta_tpos, theta_rpos, phir_e, phit_e, phir_n, phit_n, Re);

    % height of the Earth's surface above sea level where the common volume is located

    Hs = surface_altitude_cv(h, d, dt_cv)/1000.0; % in km

    Ht = (htg + h(1))/1000; % in km
    Hr = (hrg + h(end))/1000; % in km
    
    [Lbs, theta_s] = tl_troposcatter_pdr(f, dtot, Ht, Hr, theta_t, theta_r, phi_cvn, phi_cve, Gt, Gr, p, Hs);
       
end

% Calculate the final transmission loss not exceeded for p% time

Lb_pol = -5*log10(10.^(-0.2*Lbs) + 10.^(-0.2*Lbam)) + Aht + Ahr;  % eq (64)

Lb = Lb_pol(pol);


return
end