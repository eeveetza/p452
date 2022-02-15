function [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f)
%smooth_earth_heights smooth-Earth effective antenna heights according to ITU-R P.452-16
% [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta_tot, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f)
% This function derives smooth-Earth effective antenna heights according to
% Sections 4 and 5 of the Annex 2 of ITU-R P.452-16
%
% Input parameters:
% d         -   vector of terrain profile distances from Tx [0,dtot] (km)
% h         -   vector of terrain profile heigths amsl (m)
% htg, hrg  -   Tx and Rx antenna heights above ground level (m)
% ae        -   median effective Earth's radius (c.f. Eq (6a))
% f         -   frequency (GHz)
%
% Output parameters:
%
% hst, hsr     -   Tx and Rx antenna heigts of the smooth-Earth surface amsl (m)
% hstd, hsrd   -   Tx and Rx effective antenna heigts for the diffraction model (m)
% hte, hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
% hm           -   The terrain roughness parameter (m)
% dlt          -   interfering antenna horizon distance (km)
% dlr          -   Interfered-with antenna horizon distance (km)
% theta_t      -   Interfering antenna horizon elevation angle (mrad)
% theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
% theta_tot    -   Angular distance (mrad)
% pathtype     -   1 = 'los', 2 = 'transhorizon'
%
% Rev   Date        Author                          Description
% -------------------------------------------------------------------------------
% v0    15JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
% v1    15JUN16     Ivica Stevanovic, OFCOM         Modifications related to LoS path

n = length(d);

dtot = d(end);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Section 5.1.6.2

v1 = 0;
for ii = 2:n
    v1 = v1 + (d(ii)-d(ii-1))*(h(ii)+h(ii-1));  % Eq (161)
end
v2 = 0;
for ii = 2:n
    v2 = v2 + (d(ii)-d(ii-1))*( h(ii)*( 2*d(ii) + d(ii-1) ) + h(ii-1) * ( d(ii) + 2*d(ii-1) ) );  % Eq (162)
end

hst = (2*v1*dtot - v2)/dtot.^2;       % Eq (163)
hsr = (v2- v1*dtot)/dtot.^2;          % Eq (164)

% Section 5.1.6.3

HH = h - (hts*(dtot-d) + hrs*d)/dtot;  % Eq (165d)

hobs = max(HH(2:n-1));                 % Eq (165a)

alpha_obt = max( HH(2:n-1)./d(2:n-1) ); % Eq (165b)

alpha_obr = max( HH(2:n-1)./( dtot - d(2:n-1) ) ); % Eq (165c)

% Calculate provisional values for the Tx and Rx smooth surface heights

gt = alpha_obt/(alpha_obt + alpha_obr);         % Eq (166e)
gr = alpha_obr/(alpha_obt + alpha_obr);         % Eq (166f)

if hobs <= 0
    hstp = hst;                                 % Eq (166a)
    hsrp = hsr;                                 % Eq (166b)
else
    hstp = hst - hobs*gt;                       % Eq (166c)
    hsrp = hsr - hobs*gr;                       % Eq (166d)
end

% calculate the final values as required by the diffraction model

if hstp >= h(1)
    hstd = h(1);                                % Eq (167a)
else
    hstd = hstp;                                % Eq (167b)
end

if hsrp > h(end)
    hsrd = h(end);                              % Eq (167c)
else
    hsrd = hsrp;                                % Eq (167d)
    
end

% Interfering antenna horizon elevation angle and distance

ii = 2:n-1;

theta = 1000 * atan( (h(ii) - hts)./(1000 * d(ii) ) - d(ii)./(2*ae) );  % Eq (152)

%theta(theta < 0) = 0;  % condition below equation (152)

theta_t = max(theta);                           % Eq (154)


theta_td = 1000 * atan( (hrs - hts)./(1000 * dtot ) - dtot./(2*ae) );  % Eq (153)
theta_rd = 1000 * atan( (hts - hrs)./(1000 * dtot ) - dtot./(2*ae) );  % not defined in the recommendation

if theta_t > theta_td   % Eq (150): test for the trans-horizon path
    pathtype = 2; %transhorizon
else
    pathtype = 1; %los
end


kindex = find(theta == theta_t);

lt = kindex(1)+1;

dlt = d(lt);                             % Eq (155)

% Interfered-with antenna horizon elevation angle and distance

theta = 1000 * atan( (h(ii) - hrs)./(1000 * (dtot - d(ii)) ) - (dtot - d(ii))./(2*ae) );  % Eq (157)

%theta(theta < 0) = 0;

theta_r = max(theta);

kindex = find(theta == theta_r);
lr = kindex(end)+1;

dlr = dtot - d(lr);                            % Eq (158)

if pathtype == 1
   
    theta_t = theta_td;
    theta_r = theta_rd;
    
    ii = 2:n-1;
    
    lambda = 0.3/f;
    Ce = 1/ae;
    
    nu = (h(ii) + 500*Ce*d(ii).*(dtot-d(ii))- (hts*(dtot- d(ii)) + hrs *d(ii))/dtot).* ...
         sqrt(0.002*dtot./(lambda*d(ii).*(dtot-d(ii))));
    numax = max(nu);
    
    kindex = find(nu == numax);
    lt = kindex(end)+1;  
    dlt = d(lt);  
    dlr = dtot - dlt;
    kindex = find (dlr <=dtot -d(ii));
    lr = kindex(end)+1;
end

% Angular distance

theta_tot = 1e3 * dtot/ae + theta_t + theta_r;         % Eq (159)


% Section 5.1.6.4 Ducting/layer-reflection model

% Calculate the smooth-Earth heights at transmitter and receiver as
% required for the roughness factor

hst = min(hst, h(1));                           % Eq (168a)
hsr = min(hsr, h(end));                         % Eq (168b)

% Slope of the smooth-Earth surface

m = (hsr - hst)/ dtot;                          % Eq (169)

% The terminal effective heigts for the ducting/layer-reflection model

hte = htg + h(1) -   hst;                       % Eq (170)
hre = hrg + h(end) - hsr;                       

ii = lt:1:lr;

hm = max(h(ii) - (hst + m*d(ii)));

return
end