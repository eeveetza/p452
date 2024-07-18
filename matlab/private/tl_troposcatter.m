function [Lbs, theta] = tl_troposcatter(f, dt, hts, hrs, ae, the, thetat, thetar, phicvn, phicve, Gt, Gr, p, hs)
%tl_troposcatter_pdr Troposcatter basic transmission loss
%   This function computes the troposcatter basic transmission loss
%   as defined in Section 4.3 
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
%     v5    25APR23     Ivica Stevanovic, OFCOM         Aligned with rest of Recommendation 

fMHz = f*1000;  % from GHz in MHz

%% Climatic classification

% Find average annual sea-level surface refractivity N0 and radio-refractivity lapse rate dN 
% for the common volume of the link in question using the digital maps at phicve (lon),
% phicvn (lat) - as a bilinear interpolation

dN = get_interp2('DN50',phicve,phicvn);

N0 = get_interp2('N050',phicve,phicvn);


%% Calculation of tropocscatter basic transmission loss

% Step2: Calculate the scatter angle theta 

theta = 1000*the + thetat + thetar; % mrad    (145)

% Step 3: Estimate the aperture-to-median coupling loss Lc (11)

Lc = 0.07 * exp(0.055* (Gt + Gr));   % dB    (45a)

% Step 4: Estimate the average annual transmission loss associated with
% troposcatter not exceeded for p% of time (45): 

hb = 7.35;  %km  scale height set to the global mean

beta = dt/(2*ae) + thetar/1000 + (hrs-hts)/(1000*dt);  %(45e)

h0 = hts/1000 + dt*sin(beta)/(sin(theta/1000)) *(0.5* dt * sin(beta)/(ae*sin(theta/1000))+sin(thetat/1000));   %(45d)

Yp = 0.035*N0*exp(-h0/hb)*(-log10(p/50)).^(0.67);

if (p >= 50) 
    
    Yp =  -0.035*N0*exp(-h0/hb)*(-log10((100-p)/50)).^(0.67);
    
end

F = 0.18*N0*exp(-hs/hb) - 0.23*dN;   %(45b)

Lbs = F + 22.0*log10(fMHz) + 35.0*log10(theta) + 17.0*log10(dt) + Lc - Yp;    % (45)

return
end