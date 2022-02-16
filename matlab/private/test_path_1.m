function [success, fail] = test_path_1()

success = 0;
fail = 0;

disp('Performing tests for the path profile  1')
disp('as given in ITUR_452_14.')

disp(' ')


% add path to the folder where the functions are defined
s = pwd;
s=s(1:end-5);
if (~exist('longest_cont_dist.m','file'))
    addpath([s '/src'])
end
if (~exist('tl_p452.m','file'))
    addpath(s)
end


f = 0.2; % GHz
p = 0.1; % %
DN = 53;
N0 = 328;
Gt = 20;
Gr = 5;
pressure = 1013; % (hPa)
temp = 15;  % (deg C)
phi_t = 51.2;   % deg
phi_r = 50.73;  % deg
w = 0.3945;  % %
htg = 10;  % m
hrg = 10;  % m
dct = 500; %km
dcr = 500; %km

dk_t = [];
dk_r = [];
ha_t = [];
ha_r = [];

% import test profile

[d, h, zone] = test_profile(1);


% Path center latitude
phi_path = (phi_t + phi_r)/2;

% Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
zone_r = 12;
dtm = longest_cont_dist(d, zone, zone_r);

% Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
zone_r = 2;
dlm = longest_cont_dist(d, zone, zone_r);

% Compute b0
b0 = beta0(phi_path, dtm, dlm);

[ae, ab] = earth_rad_eff(DN);

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;


% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

%% test smooth earth parameters
ref = [
9.6177596154E+03
1.0900000000E+02
5.0000000000E+01
1.9300000000E+02
-6.3421180035E-01
-1.3900396738E+00
9.3089492253E+00
4.8689504250E+00
6.6222792694E+01
1.1952326473E+02
4.4582947563E+01
1.2189411666E+02
2.8000000000E+01
1.1000000000E+01
3.4500000000E+01
6.0000000000E+00
3.2827313891E+00
].';

comp = [ae, dtot, hts, hrs, theta_t, theta_r, theta, hstd, hsrd, hm, hte, hre, dlt, dlr, dtm, dlm, b0];

error = min(-log10(abs(comp-ref)));

if error > 7
    fprintf(1,'1... Smooth earth parameters, passed\n');
    success = success + 1;
else
    fprintf(1,'1... Smooth earth parameters, failed\n');
    fprintf(1,'1... Minimum significant digits agreement %f\n',error);
    fail = fail + 1;
end



%% Testing the complete propagation model

% Clutter attenuation

Lb_ref = [
135.757135679968  
139.141163124942  
135.757229717689  
130.228916886651  
125.140650939869  
129.254129344984  
133.421551585368  
137.641195863568  
141.945756530791  
146.387574001096  
151.045853314396  
156.084530033867  
162.068505261201  
176.042506912628  
182.338174670292  
197.025737048127  
    ];

ff = [
1.0000000000E-01            % 135.757135679968
1.5000000000E-01            % 139.141163124942
2.2500000000E-01            % 135.757229717689
3.3750000000E-01            % 130.228916886651
5.0625000000E-01            % 125.140650939869
7.5937500000E-01            % 129.254129344984
1.1390625000E+00            % 133.421551585368
1.7085937500E+00            % 137.641195863568
2.5628906250E+00            % 141.945756530791
3.8443359375E+00            % 146.387574001096
5.7665039063E+00            % 151.045853314396
8.6497558594E+00            % 156.084530033867
1.2974633789E+01            % 162.068505261201
1.9461950684E+01            % 176.042506912628
2.9192926025E+01            % 182.338174670292
4.3789389038E+01            % 197.025737048127
];
Aht = 0;
Ahr = 0;
pol = 1; %polarization
for ii = 1:length(ff)
Lb(ii,:) = tl_p452(ff(ii), p, d, h, zone, htg, hrg, phi_path, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp);
end


comp = Lb(:,1);
ref = Lb_ref;
error = max(abs(comp-ref));
if (error < 1.5)
    disp('2...  Frequency sweep passed')
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    success = success + 1;
else
    fprintf(1,'2...  Frequency sweep failed\n');
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail+1;
end

fprintf(1,'\n');


% plot(ff, comp,'b', ff, ref,'r--')
% legend('MATLAB','EXCEL')


%% Testing the complete propagation model scaling with p



Lb_ref = [
151.326259486794 %1
166.947740934563 %4
175.390594206107 %7
177.958970352001 %10
179.978093430805 %13
181.668993638618 %16
183.139739727119 %19
184.451368378848 % 22
185.641697122223 % 25
186.735762891508 % 28
187.751077804249 % 31
188.700617259145 % 34
189.5947770161 % 37
190.442968362355 % 40
191.255504968248 % 43
192.047496871087 %46
192.857634763625
    ];

f = 2; 

pp = [1
4
7
10
13
16
19
22
25
28
31
34
37
40
43
46
49
];
Aht = 0;
Ahr = 0;
pol = 1; %polarization
for ii = 1:length(pp)
Lb(ii,:) = tl_p452(f, pp(ii), d, h, zone, htg, hrg, phi_path, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp );
end


comp = Lb(:,1);
ref = Lb_ref;
error = max(abs(comp-ref));
if (error < 0.5)
    disp('3...  Complete propagation test scaling with p passed')
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    success = success + 1;
else
    fprintf(1,'3...  Complete propagation test scaling with p failed\n');
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail+1;
end

fprintf(1,'\n');

% figure
% xlabel('p %')
% plot(pp, comp,'b', pp, ref,'r--')
% legend('MATLAB','EXCEL')

return
end