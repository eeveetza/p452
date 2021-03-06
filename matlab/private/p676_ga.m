function [g_0, g_w] = p676_ga(f, p, rho, T)
%p676_ga Specific attenuation due to dry air and water vapour
% [g_0, g_w] = p676_ga(f, p, rho, T)
% This function computes the specific attenuation due to dry air and water vapour,
% at frequencies up to 1 000 GHz for different values of of pressure, temperature
% and humidity by means of a summation of the individual resonance lines from
% oxygen and water vapour according to ITU-R P.676-10
%
% Input parameters:
% f       -   Frequency (GHz)
% p       -   Dry air pressure (hPa)
% rho     -   Water vapor density (g/m^3)
% T       -   Temperature (K)
%
% Output parameters:
% g_o, g_w   -   specific attenuation due to dry air and water vapour
%
%
% Rev   Date        Author                          Description
% ----------------------------------------------------------------------------------------------
% v0    04FEB14     Ivica Stevanovic, OFCOM         First implementation in matlab
% v1    15DEC15     Ivica Stevanovic, OFCOM         Introduced comments and fine tuned the code
%
%

%% spectroscopic data for oxigen
%             f0        a1    a2     a3   a4     a5     a6
oxigen = [50.474214, 0.975, 9.651, 6.690, 0.0, 2.566, 6.850;
    50.987745, 2.529, 8.653, 7.170, 0.0, 2.246, 6.800;
    51.503360, 6.193, 7.709, 7.640, 0.0, 1.947, 6.729;
    52.021429, 14.320, 6.819, 8.110, 0.0, 1.667, 6.640;
    52.542418, 31.240, 5.983, 8.580, 0.0, 1.388, 6.526;
    53.066934, 64.290, 5.201, 9.060, 0.0, 1.349, 6.206;
    53.595775, 124.600, 4.474, 9.550, 0.0, 2.227, 5.085;
    54.130025, 227.300, 3.800, 9.960, 0.0, 3.170, 3.750;
    54.671180, 389.700, 3.182, 10.370, 0.0, 3.558, 2.654;
    55.221384, 627.100, 2.618, 10.890, 0.0, 2.560, 2.952;
    55.783815, 945.300, 2.109, 11.340, 0.0, -1.172, 6.135;
    56.264774, 543.400, 0.014, 17.030, 0.0, 3.525, -0.978;
    56.363399, 1331.800, 1.654, 11.890, 0.0, -2.378, 6.547;
    56.968211, 1746.600, 1.255, 12.230, 0.0, -3.545, 6.451;
    57.612486, 2120.100, 0.910, 12.620, 0.0, -5.416, 6.056;
    58.323877, 2363.700, 0.621, 12.950, 0.0, -1.932, 0.436;
    58.446588, 1442.100, 0.083, 14.910, 0.0, 6.768, -1.273;
    59.164204, 2379.900, 0.387, 13.530, 0.0, -6.561, 2.309;
    59.590983, 2090.700, 0.207, 14.080, 0.0, 6.957, -0.776;
    60.306056, 2103.400, 0.207, 14.150, 0.0, -6.395, 0.699;
    60.434778, 2438.000, 0.386, 13.390, 0.0, 6.342, -2.825;
    61.150562, 2479.500, 0.621, 12.920, 0.0, 1.014, -0.584;
    61.800158, 2275.900, 0.910, 12.630, 0.0, 5.014, -6.619;
    62.411220, 1915.400, 1.255, 12.170, 0.0, 3.029, -6.759;
    62.486253, 1503.000, 0.083, 15.130, 0.0, -4.499, 0.844;
    62.997984, 1490.200, 1.654, 11.740, 0.0, 1.856, -6.675;
    63.568526, 1078.000, 2.108, 11.340, 0.0, 0.658, -6.139;
    64.127775, 728.700, 2.617, 10.880, 0.0, -3.036, -2.895;
    64.678910, 461.300, 3.181, 10.380, 0.0, -3.968, -2.590;
    65.224078, 274.000, 3.800, 9.960, 0.0, -3.528, -3.680;
    65.764779, 153.000, 4.473, 9.550, 0.0, -2.548, -5.002;
    66.302096, 80.400, 5.200, 9.060, 0.0, -1.660, -6.091;
    66.836834, 39.800, 5.982, 8.580, 0.0, -1.680, -6.393;
    67.369601, 18.560, 6.818, 8.110, 0.0, -1.956, -6.475;
    67.900868, 8.172, 7.708, 7.640, 0.0, -2.216, -6.545;
    68.431006, 3.397, 8.652, 7.170, 0.0, -2.492, -6.600;
    68.960312, 1.334, 9.650, 6.690, 0.0, -2.773, -6.650;
    118.750334, 940.300, 0.010, 16.640, 0.0, -0.439, 0.079;
    368.498246, 67.400, 0.048, 16.400, 0.0, 0.000, 0.000;
    424.763020, 637.700, 0.044, 16.400, 0.0, 0.000, 0.000;
    487.249273, 237.400, 0.049, 16.000, 0.0, 0.000, 0.000;
    715.392902, 98.100, 0.145, 16.000, 0.0, 0.000, 0.000;
    773.839490, 572.300, 0.141, 16.200, 0.0, 0.000, 0.000;
    834.145546, 183.100, 0.145, 14.700, 0.0, 0.000, 0.000];

%% spectroscopic data for water-vapor
%            f0       b1    b2    b3   b4   b5   b6
vapor = [22.235080 0.1130 2.143 28.11 .69 4.800 1.00;
    67.803960 0.0012 8.735 28.58 .69 4.930 .82 ;
    119.995940 0.0008 8.356 29.48 .70 4.780 .79;
    183.310091 2.4200 .668 30.50 .64 5.300 .85 ;
    321.225644 0.0483 6.181 23.03 .67 4.690 .54;
    325.152919 1.4990 1.540 27.83 .68 4.850 .74;
    336.222601 0.0011 9.829 26.93 .69 4.740 .61;
    380.197372 11.5200 1.048 28.73 .54 5.380 .89;
    390.134508 0.0046 7.350 21.52 .63 4.810 .55  ;
    437.346667 0.0650 5.050 18.45 .60 4.230 .48  ;
    439.150812 0.9218 3.596 21.00 .63 4.290 .52  ;
    443.018295 0.1976 5.050 18.60 .60 4.230 .50  ;
    448.001075 10.3200 1.405 26.32 .66 4.840 .67 ;
    470.888947 0.3297 3.599 21.52 .66 4.570 .65  ;
    474.689127 1.2620 2.381 23.55 .65 4.650 .64  ;
    488.491133 0.2520 2.853 26.02 .69 5.040 .72  ;
    503.568532 0.0390 6.733 16.12 .61 3.980 .43  ;
    504.482692 0.0130 6.733 16.12 .61 4.010 .45  ;
    547.676440 9.7010 .114 26.00 .70 4.500 1.00  ;
    552.020960 14.7700 .114 26.00 .70 4.500 1.00 ;
    556.936002 487.4000 .159 32.10 .69 4.110 1.00;
    620.700807 5.0120 2.200 24.38 .71 4.680 .68  ;
    645.866155 0.0713 8.580 18.00 .60 4.000 .50  ;
    658.005280 0.3022 7.820 32.10 .69 4.140 1.00 ;
    752.033227 239.6000 .396 30.60 .68 4.090 .84 ;
    841.053973 0.0140 8.180 15.90 .33 5.760 .45  ;
    859.962313 0.1472 7.989 30.60 .68 4.090 .84  ;
    899.306675 0.0605 7.917 29.85 .68 4.530 .90  ;
    902.616173 0.0426 8.432 28.65 .70 5.100 .95  ;
    906.207325 0.1876 5.111 24.08 .70 4.700 .53  ;
    916.171582 8.3400 1.442 26.70 .70 4.780 .78  ;
    923.118427 0.0869 10.220 29.00 .70 5.000 .80 ;
    970.315022 8.9720 1.920 25.50 .64 4.940 .67  ;
    987.926764 132.1000 .258 29.85 .68 4.550 .90 ;
    1780.000000 22300.0000 .952 176.20 .50 30.500 5.00];

a1 = oxigen(:,2);
a2 = oxigen(:,3);
a3 = oxigen(:,4);
a4 = oxigen(:,5);
a5 = oxigen(:,6);
a6 = oxigen(:,7);

b1 = vapor(:,2);
b2 = vapor(:,3);
b3 = vapor(:,4);
b4 = vapor(:,5);
b5 = vapor(:,6);
b6 = vapor(:,7);



theta = 300.0/T;

e = rho * T / 216.7;        % equation (4)

%% Oxigen computation
fi = oxigen(:,1);

Si = a1 .* 1e-7 * p * theta.^3 .*exp(a2 * (1.0 - theta));       % equation (3)

df = a3 .* 1e-4 .* (p .* theta .^ (0.8-a4) + 1.1 * e * theta);  % equation (6a)

% Doppler broadening

df = sqrt( df.*df + 2.25e-6);                                   % equation (6b)

delta = (a5 + a6 * theta) * 1e-4 * (p + e) .* theta.^(0.8);     % equation (7)

Fi = f ./ fi .* (  (df - delta .* (fi - f))./( (fi - f).^2 + df.^2  ) + ...
    (df - delta .* (fi + f))./( (fi + f).^2 + df.^2  ));        % equation (5)

d = 5.6e-4 * (p + e) * theta.^(0.8);                            % equation (9)

Ndf = f * p * theta.^2 *( 6.14e-5/(d * (1 + (f/d).^2) ) + ...
    1.4e-12 * p * theta.^(1.5)/(1 + 1.9e-5 * f.^(1.5)) );       % equation (8)

% specific attenuation due to dry air (oxygen, pressure induced nitrogen
% and non-resonant Debye attenuation), equations (1-2)

if  f <= 118.750343
    g_0 = 0.182 * f * (sum(Si .* Fi) + Ndf);
else
    g_0 = 0.182 * f * (sum(Si(38:end) .* Fi(38:end)) + Ndf);
end



%% vapor computation

fi = vapor(:,1);

Si = b1 .* 1e-1 .* e .* theta.^3.5 .*exp(b2 .* (1.0 - theta));      % equation (3)

df = b3 .* 1e-4 .* (p .* theta .^ (b4) + b5 .* e .* theta.^b6);     % equation (6a)

% doppler broadening

df = 0.535 .* df + sqrt( 0.217* df.*df + 2.1316e-12 * fi.*fi/theta); % equation (6b)

delta = 0;                                                           % equation (7)

Fi = f ./ fi .* (  (df - delta .* (fi - f))./( (fi - f).^2 + df.^2  ) + ...
    (df - delta.* (fi + f))./( (fi + f).^2 + df.^2  ));              % equation (5)

% specific attenuation due to water vapour, equations (1-2)

g_w = 0.182 * f * (sum(Si .* Fi) );



return
end