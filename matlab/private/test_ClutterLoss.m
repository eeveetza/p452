function [success, fail] = test_ClutterLoss()

success = 0;
fail = 0;

disp('Performing tests for the path profile  1')
disp('as given in ITUR_452_14.')
disp('for Clutter Loss computation ')



% add path to the folder where the functions are defined

% add path to the folder where the functions are defined
s = pwd;
s=s(1:end-5);
if (~exist('longest_cont_dist.m','file'))
    addpath([s '/src'])
end
if (~exist('tl_p452.m','file'))
    addpath(s)
end

f = 2; % GHz
p = 49; % %
DN = 53;
N0 = 328;
Gt = 20;
Gr = 5;
pressure = 1013; % (hPa)
temp = 15;  % (deg C)
phi_t = 51.2;   % deg
phi_r = 50.73;  % deg
htg = 10;  % m
hrg = 10;  % m
dct = 500; %km
dcr = 500; %km

% import test profile

[d, h, zone] = test_profile(1);


pol = 1; %polarization vertical

% Industrial zone for both Rx and Tx
ha_t = [4  5 15 20 9 12 20 25 35 20];
dk_t = [0.1 0.07 0.05 0.05 0.025 0.02 0.02 0.02 0.02 0.05];
ha_r = [4  5 15 20 9 12 20 25 35 20];
dk_r = [0.1 0.07 0.05 0.05 0.025 0.02 0.02 0.02 0.02 0.05];

for ii = 1:length(ha_t)
    Lb(ii) = tl_p452(f, p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp, ha_t(ii), ha_r(ii), dk_t(ii), dk_r(ii));
end

ref = [192.852982844169 ...
       192.852982844169 ...
       209.818069521283 ... 
       225.37703794711 ...
       192.852982844169 ...
       199.190346859552 ...
       226.348103202261 ...
       229.814902341375 ...
       228.643083463793 ...
       225.37703794711];

comp = Lb;

error = max(abs(comp-ref));
if (error < 0.1)
    disp('1...  Complete propagation test with clutter losses passed')
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    success = success + 1;
else
    fprintf(1,'3...  Complete propagation test with clutter losses failed\n');
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail+1;
end

fprintf(1,'\n');

% figure
% xlabel('p %')
% plot(pp, comp,'b', pp, ref,'r--')
% legend('MATLAB','EXCEL')

return
%end