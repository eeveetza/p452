function [Lbs, theta] = tl_troposcatter_pdr(f, dt, ht, hr, thetat, thetar, phicvn, phicve,  Gt, Gr, p, hs)
%tl_troposcatter_pdr Troposcatter basic transmission loss
%   This function computes the troposcatter basic transmission loss
%   as defined in WP 3M Chairman's Report 3M/364 Annex 2 
%
%     Input parameters:
%     f       -   Frequency GHz
%     dt      -   Total distance (km)
%     ht, hr  -   Altitudes of transmitting antenna and receiving antennas in km  
%     thetat  -   Tx horizon elevation angle relative to the local horizontal (mrad)
%     thetar  -   Rx horizon elevation angle relative to the local horizontal (mrad)
%     phicvn  -   Troposcatter common volume latitude (deg)
%     phicve  -   Troposcatter common volume longitude (deg)
%     Gt, Gr  -   Gain of transmitting and receiving antenna in the azimuthal direction
%                 of the path towards the other antenna and at the elevation angle
%                 above the local horizontal of the other antenna in the case of a LoS
%                 path, otherwise of the antenna's radio horizon, for median effective
%                 Earth radius.
%     p       -   Percentage of average year for which predicted basic loss
%                 is not exceeded (%)
%     hs      -   height of the earth's surface above sea level (km) 
%
%     Output parameters:
%     Lbs    -   Troposcatter basic transmission loss (dB)
%     theta  -   Scatter angle (mrad)
%
%     Example:
%     [Lbs, theta] = tl_troposcatter_pdr(f, dt, ht, hr, thetat, thetar, phicvn, phicve,  Gt, Gr, p, hs)

%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    18JUL16     Ivica Stevanovic, OFCOM         Initial version
%     v1    13JUN17     Ivica Stevanovic, OFCOM         replaced load function calls to increase computational speed
%     v2    11AUG17     Ivica Stevanovic, OFCOM         introduced a correction in TropoClim vector indices to cover the case 
%                                                       when the point is equally spaced from the two closest grid points
%     v3    28OCT19     Ivica Stevanovic, OFCOM         Introduced proposal from China as given in Annex 8 of the Chairmans Report 3M/343-E from Montreal meeting 2018
%     v4    15FEB23     Ivica Stevanovic, OFCOM         Introduced the updates from PDR (see WP 3M Chairman's report 2022 3M/364 Annex 2 )


% Attachment E: Troposcatter

fMHz = f*1000;  % from GHz in MHz

%% E.2 Climatic classification

DN50 = DigitalMaps_DN50();

N050 = DigitalMaps_N050();

latcnt = 90:-1.5:-90;               %Table 2.4.1
loncnt = 0:1.5:360;                 %Table 2.4.1  
[LON,LAT] = meshgrid(loncnt, latcnt);

% Map phicve (-180, 180) to loncnt (0,360);

phicve1 = phicve;
if phicve1 < 0
    phicve1 = phicve + 360;
end

% Find average annual sea-level surface refractivity N0 and radio-refractivity lapse rate dN 
% for the common volume of the link in question using the digital maps at phicve (lon),
% phicvn (lat) - as a bilinear interpolation

dN = interp2(LON,LAT,DN50,phicve1,phicvn);

N0 = interp2(LON,LAT,N050,phicve1,phicvn);

clear DN50 N050


%% E.3 Calculation of tropocscatter basic transmission loss

% Step2: Calculate the scatter angle theta 
k = 4.0/3.0;
a = 6371;
the = dt * 1000 /(k * a);   % mrad   (E.2)

theta = the + thetat + thetar; % mrad    (E.1)

% Step 3: Estimate the aperture-to-median coupling loss Lc (11)

Lc = 0.07 * exp(0.055* (Gt + Gr));   % dB    (E.3)

% Step 4: Estimate the average annual transmission loss associated with
% troposcatter not exceeded for p% of time (E.4)-(E.8): 

hb = 7.35;  %km  scale height set to the global mean

beta = dt/(2*k*a) + thetar/1000 + (hr-ht)/dt;  %(E.8)

h0 = ht + dt*sin(beta)/(sin(theta/1000)) *(0.5* dt * sin(beta)/(k*a*sin(theta/1000))+sin(thetat/1000));   %(E.7)

Yp = 0.035*N0*exp(-h0/hb)*(-log10(p/50)).^(0.67);

if (p >= 50) 
    
    Yp =  -0.035*N0*exp(-h0/hb)*(-log10((100-p)/50)).^(0.67);
    
end

F = 0.18*N0*exp(-hs/hb) - 0.23*dN;   %(E.5)

Lbs = F + 22.0*log10(fMHz) + 35.0*log10(theta) + 17.0*log10(dt) + Lc - Yp;    % (E.4)

return
end