function [succ, fail] = test_path_2()

disp('Performing tests for the path profile  2')
disp('as given in ITUR_452_15_REV4.')
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

succ = 0;
fail = 0;

f = 1.8; % GHz
p = 10; % %
DN = 50;
N0 = 301;
Gt = 12.3796;
Gr = 22.0409;
pressure = 1013; % (hPa)
temp = 15;  % (deg C)
phi_t = 40.5;   % deg
phi_r = 40;  % deg
w = 0.3945;  % %
htg = 10;  % m
hrg = 10;  % m
dct = 500; %km
dcr = 500; %km
thshld = 0.8;
% import test profile

[d, h, zone] = test_profile(2);


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



%% test line-of-sight propagation

ff = [1.0000000000E-01
1.5000000000E-01
2.2500000000E-01
3.3750000000E-01
5.0625000000E-01
7.5937500000E-01
1.1390625000E+00
1.7085937500E+00
2.5628906250E+00
3.8443359375E+00
5.7665039063E+00
8.6497558594E+00
1.2974633789E+01
1.9461950684E+01
2.9192926025E+01
4.3789389038E+01
].';

for ii = 1:length(ff)
    
[Lbfsg(ii), Lb0p(ii), Lb0b(ii)] = pl_los(dtot, ff(ii), p, b0, omega, temp, pressure, dlt, dlr);

end

ref = [1.0940899743E+02
1.1294774579E+02
1.1650370095E+02
1.2008691031E+02
1.2369942906E+02
1.2732301193E+02
1.3093023162E+02
1.3451125308E+02
1.3807708122E+02
1.4164957136E+02
1.4526463082E+02
1.4900478082E+02
1.5317811739E+02
1.6177617665E+02
1.6589663737E+02
1.7504054471E+02
].';

comp = Lbfsg;

error = max((abs(comp-ref)));

if (error < thshld)
    disp('2.1..  Line-of-sight frequency sweep test (Lbfsg) passed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    succ = succ+1;
else
    fprintf(1,'2.1..  Line-of-sight frequency sweep test (Lbfsg) failed\n');
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    fail = fail+ 1;
end
 

ref = [108.2329949
111.7717432
115.3276984
118.9109077
122.5234265
126.1470094
129.754229
133.3352505
136.9010786
140.4735688
144.0886282
147.8287782
152.0021148
160.6001741
164.7206348
173.8645421
].';

comp = Lb0p;

error = max((abs(comp-ref)));

if (error < thshld)
    disp('2.2..  Line-of-sight frequency sweep test (Lb0p) passed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    succ = succ + 1;
else
    fprintf(1,'2.2..  Line-of-sight frequency sweep test (Lb0p) failed\n');
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    fail = fail + 1;
end


ref = [1.0723849158E+02
1.1077723994E+02
1.1433319510E+02
1.1791640446E+02
1.2152892321E+02
1.2515250608E+02
1.2875972577E+02
1.3234074723E+02
1.3590657537E+02
1.3947906551E+02
1.4309412497E+02
1.4683427497E+02
1.5100761154E+02
1.5960567080E+02
1.6372613152E+02
1.7287003886E+02
].';

comp = Lb0b;

error = max((abs(comp-ref)));

if (error < thshld)
    disp('2.3..  Line-of-sight frequency sweep test (Lb0b) passed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    succ = succ + 1;
else
    fprintf(1,'2.3..  Line-of-sight frequency sweep test (Lb0b) failed\n');
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    fail = fail + 1;
end


fprintf(1,'\n');
 
%% Test troposcatterer transmission loss
for ii = 1:length(ff)
Lbs(ii) = tl_tropo(dtot, theta, ff(ii), p, temp, pressure, N0, Gt, Gr );
end
comp = Lbs.';
ref = [1.5929883240E+02
1.6478598472E+02
1.7013526311E+02
1.7535666489E+02
1.8045213249E+02
1.8540316673E+02
1.9018176773E+02
1.9477680524E+02
1.9919633333E+02
2.0345536835E+02
2.0757327763E+02
2.1158762074E+02
2.1562970174E+02
2.2134185316E+02
2.2528794965E+02
2.3332411115E+02
];

error = max((abs(comp-ref)));



if (error < thshld)
    disp('4...  Tropospheric frequency sweep test passed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    succ = succ + 1;
else
    fprintf(1, '4...  Tropospheric frequency sweep test failed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    fail = fail + 1;
end


fprintf(1,'\n');
%% Test anomalous transmission loss
for ii = 1: length(ff)
Lba(ii) = tl_anomalous(dtot, dlt, dlr, dct, dcr, dlm, hts, hrs, hte, hre, hm, theta_t, theta_r, ff(ii), p, temp, pressure, w, ae, b0);
end

comp = Lba.';
ref = [184.0635635
183.6901287
181.5069767
177.623325
173.818353
179.8708527
186.0904005
192.4809533
199.0707418
205.9027522
213.0384176
220.5907073
228.9036759
242.0109963
251.0594292
265.6069331
];

error = max((abs(comp-ref)));

if (error < 2)
    disp('5...  Anomalous frequency sweep test passed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    succ = succ + 1;
else
    fprintf(1, '5...  Anomalous frequency sweep test failed')
    fprintf(1,'     (Maximum distance from reference: %f dB)\n', error);
    fail = fail + 1;
end

fprintf(1,'\n');

%% Testing the complete propagation model

% Clutter attenuation

Lb_ref = [
    156.097581791027
    159.479671293708
    162.905190981865
    166.568770512145
    170.511694572475
    174.762717586798
    179.283540982975
    183.973934469392
    188.804623137943
    193.905732182431
    199.438463405793
    205.438453008016
    212.020288608104
    220.828620828485
    225.18672608976
    233.314425507872
    ];

ff = [1.0000000000E-01      % 156.097581791027
1.5000000000E-01            % 159.479671293708
2.2500000000E-01            % 162.905190981865
3.3750000000E-01            % 166.568770512145
5.0625000000E-01            % 170.511694572475
7.5937500000E-01            % 174.762717586798
1.1390625000E+00            % 179.283540982975
1.7085937500E+00            % 183.973934469392
2.5628906250E+00            % 188.804623137943
3.8443359375E+00            % 193.905732182431
5.7665039063E+00            % 199.438463405793
8.6497558594E+00            % 205.438453008016
1.2974633789E+01            % 212.020288608104
1.9461950684E+01            % 220.828620828485
2.9192926025E+01            % 225.18672608976
4.3789389038E+01            % 233.314425507872
];

pol = 1; %polarization
for ii = 1:length(ff)
Lb(ii,:) = tl_p452(ff(ii), p, d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp );
end


comp = Lb(:,1);
ref = Lb_ref;
error = max(abs(comp-ref));
if (error < thshld)
    disp('6.1...  Complete propagation test passed')
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'6.1...  Complete propagation test failed\n');
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

fprintf(1,'\n');
% 
% 
% plot(ff, comp,'b', ff, ref,'r--')
% legend('MATLAB','EXCEL')

%% Testing the complete propagation model scaling with p

% Clutter attenuation

Lb_ref = [
    163.0933482
178.1843632
184.4675066
185.8333773
186.9079496
187.8103368
188.5990665
189.3075163
189.9566814
190.5607394
191.1298072
191.6714433
192.1915287
192.6948222
193.1853415
193.6666828
194.142628
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

pol = 1; %polarization
for ii = 1:length(pp)
Lb(ii,:) = tl_p452(f, pp(ii), d, h, zone, htg, hrg, phi_t, phi_r, Gt, Gr, pol, dct, dcr, DN, N0, pressure, temp );
end


comp = Lb(:,1);
ref = Lb_ref;
error = max(abs(comp-ref));
if (error < thshld)
    disp('6.2...  Complete propagation test scaling with p passed')
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'6.2...  Complete propagation test scaling with p failed\n');
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

fprintf(1,'\n');
return
end
% figure
% xlabel('p %')
% plot(pp, comp,'b', pp, ref,'r--')
% legend('MATLAB','EXCEL')