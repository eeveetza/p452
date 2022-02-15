function Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega )
%dl_delta_bull Complete 'delta-Bullington' diffraction loss model P.452-16
%   function Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega )
%
%   This function computes the complete 'delta-Bullington' diffraction loss
%   as defined in ITU-R P.452-16 (Section 4.2.3)
%
%     Input parameters:
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%     above mean sea level. Both vectors contain n+1 profile points
%     hts     -   transmitter antenna height in meters above sea level (i=0)
%     hrs     -   receiver antenna height in meters above sea level (i=n)
%     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.1.6.3
%     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.1.6.3
%     ap      -   the effective Earth radius in kilometers
%     f       -   frequency expressed in GHz
%     omega   -   the fraction of the path over sea
%
%     Output parameters:
%     Ld     -   diffraction loss for the general patha according to
%                Section 4.2.3 of ITU-R P.452-16. 
%                Ld(1) is for the horizontal polarization 
%                Ld(2) is for the vertical polarization
%
%     Example:
%     Ld = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega )
%       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version


%% Body of function

% Use the method in 4.2.1 for the actual terrain profile and antenna
% heights. Set the resulting Bullington diffraction loss for the actual
% path to Lbulla

Lbulla = dl_bull(d, h, hts, hrs, ap, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ap, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ap, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)

return
end