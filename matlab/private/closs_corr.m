function  [dc, hc, zonec,htgc, hrgc, Aht, Ahr] = closs_corr(f, d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r)
%closs clutter loss correction according to P.452-16
%   function [dc,hc,zonec,htgc,htrc, Aht, Ahr] = closs_corr(d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r)
%
%   This function computes the height-gain correction as defined in ITU-R P.452-16 (Section 4.5.4)
%
%     Input parameters:
%     f       -   Frequency (GHz)
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%                 above mean sea level. Both vectors contain n+1 profile points
%     zone    -   Zone type: Coastal land (1), Inland (2) or Sea (3)
%     htg     -   Tx Antenna center heigth above ground level (m)
%     hrg     -   Rx Antenna center heigth above ground level (m)
%     ha_t    -   Nominal clutter height at the transmitting end (m, agl)
%     ha_r    -   Nominal clutter height at the receiving end (m, agl)
%     dk_t    -   distance from nominal clutter point to the Tx antenna (km)
%     dk_r    -   distance from nominal clutter point to the Rx antenna (km)
%
%     Output parameters:
%     dc      -   vector of distances in the height-gain model
%     hc      -   vector of heights in the height-gain model
%     zonec   -   Zone type: Coastal land (1), Inland (2) or Sea (3)in the height-gain model
%     htgc    -   Tx Antenna center heigth above ground level (m) in the height-gain model
%     hrgc    -   Rx Antenna center heigth above ground level (m) in the height-gain model
%     Aht     -   Additional losses to account for clutter shielding the
%     Ahr         transmitter and receiver. These should be set to zero if there is no
%                 such shielding
%
%     Example:
%     [dc,hc,zonec,htgc,hrgc, Aht, Ahr] = closs_corr(d, h, zone, htg, hrg, ha_t, ha_r, dk_t, dk_r)
%
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    09MAR16     Ivica Stevanovic, OFCOM         Initial version
%     v1    24NOV16     Ivica Stevanovic, OFCOM         introduced clutter nominal distance check


index1 = 1;
index2 = length(d);

htgc = htg;
hrgc = hrg;

Aht = 0;
Ahr = 0;

ha = ha_t;
dk = dk_t;

if ha > htg
    
    Ffc = 0.25+0.375*(1+tanh(7.5*(f-0.5)));  % (57a)
    
    Aht = 10.25*Ffc*exp(-dk)*( 1- tanh(6*(htg/ha-0.625)) )-0.33; % (57)
    
    flagAht = 1;
    
    kk = find(d>=dk);
    
    if (~isempty(kk))
        index1 = kk(1);
    else
        index1 = length(d);
    end
    
    htgc = ha_t;

    
end

ha = ha_r;
dk = dk_r;

if ha > hrg
    
    Ffc = 0.25+0.375*(1+tanh(7.5*(f-0.5)));  % (57a)
    
    Ahr = 10.25*Ffc*exp(-dk)*( 1- tanh(6*(hrg/ha-0.625)) )-0.33;  % (57)
    
    flagAhr = 1;
    
    kk = find(d <= d(end)-dk);
    if(~isempty(kk))
        index2 = kk(end);
    else
        index2 = 1;
    end
    
    hrgc = ha_r;
end

% Modify the path

if (index2-index1 < 4) % at least three points between the clutter at Tx and Rx sides
    error('tl_p452: closs_corr: the sum of clutter nominal distances is larger than the path length.');
end
    
dc = d(index1:index2)-d(index1);
hc = h(index1:index2);
zonec = zone(index1:index2);

return
end
 
