function Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, varargin)
%tl_p452 basic transmission loss according to P.452-16
%   Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, ha_t, ha_r, dk_t, dk_r )
%
%   This is the MAIN function that computes the basic transmission loss not exceeded for p% of time
%   as defined in ITU-R P.452-16 (Section 4.6). Other functions called from
%   this function are in ./src/ subfolder.
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
%     phi_t   -   Latitude of Tx station (degrees)
%     phi_r   -   Latitude of Rx station (degrees)
%     Gt, Gr  -   Antenna gain in the direction of the horizon along the
%                 great-circle interference path (dBi)
%     pol     -   polarization of the signal (1) horizontal, (2) vertical
%     dct     -   Distance over land from the transmit and receive
%     dcr         antennas to the coast along the great-circle interference path (km).
%                 Set to zero for a terminal on a ship or sea platform
%     DN      -   The average radio-refractive index lapse-rate through the
%                 lowest 1 km of the atmosphere (it is a positive quantity in this
%                 procedure) (N-units/km)
%     N0      -   The sea-level surface refractivity, is used only by the
%                 troposcatter model as a measure of location variability of the
%                 troposcatter mechanism. The correct values of DN and N0 are given by
%                 the path-centre values as derived from the appropriate
%                 maps (N-units)
%     press   -   Dry air pressure (hPa)
%     temp    -   Air temperature (degrees C)
%     ha_t    -   Clutter nominal height (m) at the Tx side
%     ha_r    -   Clutter nominal height (m) at the Rx side
%     dk_t    -   Clutter nominal distance (km) at the Tx side
%     dk_r    -   Clutter nominal distance (km) at the Rx side
%
%     Output parameters:
%     Lb     -   basic  transmission loss according to P.452-16
%
%     Example:
%     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp)
%     Lb = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, press, temp, ha_t, ha_r, dk_t, dk_r)

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

%   


% MATLAB Version 8.3.0.532 (R2014a) used in development of this code
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
% IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
%
% THE AUTHOR(S) AND OFCOM (CH) DO NOT PROVIDE ANY SUPPORT FOR THIS SOFTWARE
%
% This function calls other functions that are placed in the ./src folder
% Test functions to verify/validate the current implementation are placed  in ./test folder

s = pwd;
if ~exist('p676d11_ga.m','file')
    addpath([s '/src/'])
end

% Read the input arguments 

if nargin > 22    warning(strcat('tl_p452: Too many input arguments; The function requires at most 22',...
        'input arguments. Additional values ignored. Input values may be wrongly assigned.'));
end

if nargin <18 
    error('tl_p452: function requires at least 18 input parameters.');
end

ha_t = [];
ha_r = [];
dk_t = [];
dk_r = [];

narg = 19;


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

% Compute the path profile parameters
% Path center latitude
phi_path = (phi_t + phi_r)/2;

% Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
zone_r = 12;
dtm = longest_cont_dist(d, zone, zone_r);

% Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
zone_r = 2;
dlm = longest_cont_dist(d, zone, zone_r);

% Compute b0
b0 = beta0(phi_path, dtm, dlm);

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

% Wavelength in meters

lambda = 0.3/f;


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

[Lbfsg, Lb0p, Lb0b] = pl_los(dtot, f, p, b0, omega, temp, press, dlt, dlr);


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

% Calculate the basic transmission loss due to troposcatter not exceeded
% for any time percantage p 

Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr );

% Calculate the final transmission loss not exceeded for p% time

Lb_pol = -5*log10(10.^(-0.2*Lbs) + 10.^(-0.2*Lbam)) + Aht + Ahr;  % eq (64)

Lb = Lb_pol(pol);


return
end