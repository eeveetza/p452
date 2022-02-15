function [Lbfsg, Lb0p, Lb0b] = pl_los(d, f, p, b0, w, temp, press, dlt, dlr)
%pl_los Line-of-sight transmission loss according to ITU-R P.452-16
%     This function computes line-of-sight transmission loss (including short-term effects)
%     as defined in ITU-R P.452-16.
%
%     Input parameters:
%     d       -   Great-circle path distance (km)
%     f       -   Frequency (GHz)
%     p       -   Required time percentage(s) for which the calculated basic
%                 transmission loss is not exceeded (%)
%     b0      -   Point incidence of anomalous propagation for the path
%                 central location (%)
%     w       -   Fraction of the total path over water (%)
%     temp    -   Temperature (degrees C)
%     press   -   Dry air pressure (hPa)
%     dlt     -   For a transhorizon path, distance from the transmit antenna to
%                 its horizon (km). For a LoS path, each is set to the distance
%                 from the terminal to the profile point identified as the Bullington
%                 point in the diffraction method for 50% time
%     dlr     -   For a transhorizon path, distance from the receive antenna to
%                 its horizon (km). The same note as for dlt applies here.
%
%     Output parameters:
%     Lbfsg   -   Basic transmission loss due to free-space propagation and
%                 attenuation by atmospheric gases
%     Lb0p    -   Basic transmission loss not exceeded for time percentage, p%, due to LoS propagation
%     Lb0b    -   Basic transmission loss not exceedd for time percentage, b0%, due to LoS propagation
%
%     Example:
%     [Lbfsg, Lb0p, Lb0b] = pl_los(d, f, p, b0, w, temp, press, dlt, dlr)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    04FEB14     Ivica Stevanovic, OFCOM         First implementation in matlab

T = temp + 273.15;

% water vapor density
rho = 7.5 + 2.5 * w;  % (9a)

% compute specific attenuation due to dry air and water vapor:
[g_0, g_w] = p676d11_ga(f, press, rho, T);

Ag = (g_0 + g_w) * d;  %(9)

% Basic transmission loss due to free-space propagation and attenuation
% by atmospheric gases
Lbfsg = 92.5 + 20.0*log10(f) + 20.0*log10(d) + Ag;  % (8)

% Corrections for multipath and focusing effects at p and b0
Esp = 2.6 * (1 - exp(-0.1 * (dlt + dlr) ) ) * log10(p/50);   %(10a)
Esb = 2.6 * (1 - exp(-0.1 * (dlt + dlr) ) ) * log10(b0/50);  %(10b)

% Basic transmission loss not exceeded for time percentage p% due to
% LoS propagation
Lb0p = Lbfsg + Esp;    %(11)

% Basic transmission loss not exceeded for time percentage b0% due to
% LoS propagation
Lb0b = Lbfsg + Esb;    %(12)

return
end