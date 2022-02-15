function Lbull = dl_bull(d, h, hts, hrs, ap, f)
%dl_bull Bullington part of the diffraction loss according to P.452-16
%   This function computes the Bullington part of the diffraction loss
%   as defined in ITU-R P.452-16 in 4.2.1
%
%     Input parameters:
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%     above mean sea level. Both vectors contain n+1 profile points
%     hts     -   transmitter antenna height in meters above sea level (i=0)
%     hrs     -   receiver antenna height in meters above sea level (i=n)
%     ap      -   the effective earth radius in kilometers
%     f       -   frequency expressed in GHz
%
%     Output parameters:
%     Lbull   -   Bullington diffraction loss for a given path
%
%     Example:
%     Lbull = dl_bull(d, h, hts, hrs, ap, f)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation in matlab


%% Body of function

% Effective Earth curvature Ce (km^-1)

Ce = 1/ap;

% Wavelength in meters

lambda = 0.3/f;

% Complete path length

dtot = d(end)-d(1);

% Find the intermediate profile point with the highest slope of the line
% from the transmitter to the point

di = d(2:end-1);
hi = h(2:end-1);

Stim = max((hi + 500*Ce*di.*(dtot - di) - hts)./di );           % Eq (14)

% Calculate the slope of the line from transmitter to receiver assuming a
% LoS path

Str = (hrs - hts)/dtot;                                         % Eq (15)

if Stim < Str % Case 1, Path is LoS
    
    % Find the intermediate profile point with the highest diffraction
    % parameter nu:
    
    numax = max (...
                  ( hi + 500*Ce*di.*(dtot - di) - ( hts*(dtot - di) + hrs*di)/dtot ) .* ...
                   sqrt(0.002*dtot./(lambda*di.*(dtot-di))) ...
                 );
              
    Luc = 0;
    if numax > -0.78
        Luc = 6.9 + 20*log10(sqrt((numax-0.1).^2+1) + numax - 0.1);   % Eq (13), (17)
    end
else
    
    % Path is transhorizon
    
    % Find the intermediate profile point with the highest slope of the
    % line from the receiver to the point
    
    Srim = max((hi + 500*Ce*di.*(dtot-di)-hrs)./(dtot-di));     % Eq (18)
    
    % Calculate the distance of the Bullington point from the transmitter:
    
    dbp = (hrs - hts + Srim*dtot)/(Stim + Srim);                % Eq (19)
    
    % Calculate the diffraction parameter, nub, for the Bullington point
    
    nub =  ( hts + Stim*dbp - ( hts*(dtot - dbp) + hrs*dbp)/dtot ) * ...
                   sqrt(0.002*dtot/(lambda*dbp*(dtot-dbp)));    % Eq (20)
    
    % The knife-edge loss for the Bullington point is given by
              
    Luc = 0;
    if nub > -0.78
        Luc = 6.9 + 20*log10(sqrt((nub-0.1).^2+1) + nub - 0.1);   % Eq (13), (21)
    end
    
end

% For Luc calculated using either (17) or (21), Bullington diffraction loss
% for the path is given by

Lbull = Luc + (1 - exp(-Luc/6.0))*(10+0.02*dtot);         % Eq (22)
return
end
