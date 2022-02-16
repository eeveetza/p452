function Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
%tl_tropo Basic transmission loss due to troposcatterer to P.452-16
%   Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
%
%   This function computes the basic transmission loss due to troposcatterer 
%   not exceeded for p% of time
%   as defined in ITU-R P.452-16 (Section 4.3)
%
%     Input parameters:
%     dtot    -   Great-circle path distance (km)
%     theta   -   Path angular distance (mrad)
%     f       -   frequency expressed in GHz
%     p       -   percentage of time
%     temp    -   Temperature (deg C)
%     press   -   Dry air pressure (hPa)
%     N0      -   path centre sea-level surface refractivity derived from Fig. 6
%     Gt,Gr   -   Antenna gain in the direction of the horizon along the
%                 great-circle interference path (dBi)
%
%     Output parameters:
%     Lbs    -   the basic transmission loss due to troposcatterer 
%                not exceeded for p% of time
%
%     Example:
%     Lbs = tl_tropo(dtot, theta, f, p, temp, press, N0, Gt, Gr )
%       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version


%% Body of function

T = temp + 273.15; 

% Frequency dependent loss

Lf = 25*log10(f) - 2.5*(log10(f/2)).^2;    % eq (45a)

% aperture to medium coupling loss (dB)

Lc = 0.051*exp(0.055*(Gt+Gr));             % eq (45b)

% gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
% whole path length

rho = 3;

% compute specific attenuation due to dry air and water vapor:
[g_0, g_w] = p676d11_ga(f, press, rho, T);

Ag = (g_0 + g_w) * dtot;  %(9)

% the basic transmission loss due to troposcatter not exceeded for any time
% percentage p, below 50% is given

Lbs = 190 + Lf + 20*log10(dtot) + 0.573*theta - 0.15*N0 + Lc + Ag - 10.1*(-log10(p/50)).^(0.7);

return
end
