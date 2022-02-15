function Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, ae, b0)
%tl_anomalous Basic transmission loss due to anomalous propagation according to P.452-16
%   Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, ae, b0)
%
%   This function computes the basic transmission loss occuring during
%   periods of anomalous propagation (ducting and layer reflection)
%   as defined in ITU-R P.452-16 (Section 4.4)
%
%     Input parameters:
%     dtot         -   Great-circle path distance (km)
%     dlt          -   interfering antenna horizon distance (km)
%     dlr          -   Interfered-with antenna horizon distance (km)
%     dct, dcr     -   Distance over land from the transmit and receive
%                      antennas tothe coast along the great-circle interference path (km).
%                      Set to zero for a terminal on a ship or sea platform
%     dlm          -   the longest continuous inland section of the great-circle path (km)
%     hts, hrs     -   Tx and Rx antenna heights aobe mean sea level amsl (m)
%     hte, hre     -   Tx and Rx terminal effective heights for the ducting/layer reflection model (m)
%     hm           -   The terrain roughness parameter (m)
%     theta_t      -   Interfering antenna horizon elevation angle (mrad)
%     theta_r      -   Interfered-with antenna horizon elevation angle (mrad)
%     f            -   frequency expressed in GHz
%     p            -   percentage of time
%     temp         -   Temperature (deg C)
%     press        -   Dry air pressure (hPa)
%     omega        -   fraction of the total path over water
%     ae           -   the median effective Earth radius (km)
%     b0           -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
%
%     Output parameters:
%     Lba    -   the basic transmission loss due to anomalous propagation
%               (ducting and layer reflection)
%
%     Example:
%     Lba = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, f, p, temp, press, omega, b0)
%       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version
%     v1    10NOV16     Ivica Stevanovic, OFCOM         added a line after eq. (54) to avoid division by zero
%     v2    13FEB17     Ivica Stevanovic, OFCOM         included lower limit for alpha and upper limit for mu2


%% Body of function

% empirical correction to account for the increasing attenuation with
% wavelength inducted propagation (47a)

Alf = 0;

if f < 0.5
    Alf = 45.375 - 137.0*f + 92.5*f*f;
end

% site-shielding diffraction losses for the interfering and interfered-with
% stations (48)

theta_t1 = theta_t - 0.1*dlt;    % eq (48a)
theta_r1 = theta_r - 0.1*dlr;

Ast = 0;
Asr = 0;

if theta_t1 > 0
    Ast = 20*log10(1 + 0.361*theta_t1*sqrt(f*dlt)) + 0.264*theta_t1*f.^(1/3);
end

if theta_r1 > 0
    Asr = 20*log10(1 + 0.361*theta_r1*sqrt(f*dlr)) + 0.264*theta_r1*f.^(1/3);
end

% over-sea surface duct coupling correction for the interfering and
% interfered-with stations (49) and (49a)

Act = 0;
Acr = 0;

if dct <= 5
    if dct <= dlt
        if omega >= 0.75
            Act = -3*exp(-0.25*dct*dct)*(1+ tanh( 0.07*(50-hts) ));
        end
    end
end

if dcr <= 5
    if dcr <= dlr
        if omega >= 0.75
            Acr = -3*exp(-0.25*dcr*dcr)*(1+ tanh( 0.07*(50-hrs) ));
        end
    end
end

% specific attenuation (51)

gamma_d = 5e-5 * ae * f.^(1/3);

% angular distance (corrected where appropriate) (52-52a)

theta_t1 = theta_t;
theta_r1 = theta_r;

if theta_t> 0.1*dlt
    theta_t1 = 0.1*dlt;
end

if theta_r > 0.1*dlr
    theta_r1 = 0.1*dlr;
end

theta1 = 1e3*dtot/ae + theta_t1 + theta_r1;   

dI = min(dtot - dlt - dlr, 40);   % eq (56a)

mu3 = 1;

if hm > 10
    
    mu3 = exp( -4.6e-5 * (hm-10)*(43+6*dI) );  % eq (56)

    
end

tau = 1- exp(-(4.12e-4*dlm.^2.41));       % eq (3a)

epsilon = 3.5;

alpha = -0.6 - epsilon*1e-9*dtot.^(3.1)*tau;   % eq (55a)

if alpha < -3.4
    alpha = -3.4;
end
% correction for path geometry:

mu2 = ( 500/ae * dtot.^2/( sqrt(hte) + sqrt(hre) ).^2 ).^alpha;

if mu2 > 1
    mu2 = 1;
end

beta = b0 * mu2 * mu3;      % eq (54)

%beta = max(beta, eps);      % to avoid division by zero

Gamma = 1.076/(2.0058-log10(beta)).^1.012 * ...
        exp( -( 9.51 - 4.8*log10(beta) + 0.198*(log10(beta)).^2)*1e-6*dtot^(1.13) );

% time percentage variablity (cumulative distribution):

Ap = -12 + (1.2 + 3.7e-3*dtot)*log10(p/beta) + 12 * (p/beta).^Gamma;  % eq (53)

% time percentage and angular-distance dependent losses within the
% anomalous propagation mechanism

Adp = gamma_d*theta1 + Ap;   % eq (50)

% gaseous absorbtion derived from equation (9) using rho = 3 g/m^3 for the
% whole path length

% water vapor density

rho = 7.5 + 2.5*omega;

T = temp + 273.15;

% compute specific attenuation due to dry air and water vapor:
[g_0, g_w] = p676d11_ga(f, press, rho, T);

Ag = (g_0 + g_w) * dtot;  %(9)

% total of fixed coupling losses (except for local clutter losses) between
% the antennas and the anomalous propagation structure within the
% atmosphere (47)

Af = 102.45 + 20*log10(f) + 20*log10(dlt + dlr) + Alf + Ast + Asr + Act + Acr;

% total basic transmission loss occuring during periods of anomalaous
% propagation 

Lba = Af + Adp + Ag;

return
end