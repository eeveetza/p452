function [succ, fail] = test_deltaBullington()

disp('Performing tests ')
disp('as given in ITU-R study group validation examples')
disp('for the delta Bullington diffraction prediction method')
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

f = 1; % GHz
htg = 30;  % m
hrg = 30;  % m

p = 10; % %
DN = 39; % supposed this value, it was not given in the validation examples
N0 = 328; %
Gt = 20;
Gr = 5;
pressure = 1013; % (hPa)
temp = 15;  % (deg C)
phi_t = 40.5;   % deg
phi_r = 40;  % deg
w = 0.0;  % %

dct = 500; %km
dcr = 500; %km

% import test profile

[d, h, zone] = test_profile(3);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [32.84284766	13.74840484	15.3485837	34.44302652];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('1...  Path 1, htg = %.1f, hrg = %.1f f = %.1f passed\n',htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'1...  Path 1, htg = %f, hrg = %f f = %.1f failed\n',htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end


%% Path 1, htg = 50, hrg = 10, f = 2.5
step = 2;
path = 1;
htg = 50;
hrg = 10;
f = 2.5;
[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [37.229036	18.5970601	23.6973069	42.32928281];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%% 
step = 3;
path = 1;
htg = 20;
hrg = 20;
f = 0.6;
[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [31.22595813	15.79373777	20.91384864	36.34606901];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%% 
step = 4;
path = 1;
htg = 40;
hrg = 50;
f = 0.2;
[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [24.87792613	11.50981857	14.32933641	27.69744396

];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end



%% 
step = 5;
path = 1;
htg = 70;
hrg = 5;
f = 0.15;

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [24.78223928	14.93033347	33.13883693	42.99074274
];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end



%%
% import test profile

step = 6;
path = 2;
htg = 30;
hrg = 30;
f = 1;

[d, h, zone] = test_profile(4);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [36.20419305	35.82816509	69.45434151	69.83036947
];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 7;
path = 2;
htg = 50;
hrg = 10;
f = 2.5;

[d, h, zone] = test_profile(4);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%

ref = [44.81240928	39.66871934	85.53863749	90.68232743

];

comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 8;
path = 2;
htg = 20;
hrg = 20;
f = 0.6;
ref = [34.43917939	34.29503068	66.31647651	66.46062522
];

[d, h, zone] = test_profile(4);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end


%%
% import test profile

step = 9;
path = 2;
htg = 40;
hrg = 50;
f = 0.2;
ref = [28.37368149	27.74009568	45.7547348	46.38832061
];

[d, h, zone] = test_profile(4);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 10;
path = 2;
htg = 70;
hrg = 5;
f = 0.15;
ref = [34.61746032	26.58559997	43.35361839	51.38547874

];

[d, h, zone] = test_profile(4);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 11;
path = 3;
htg = 30;
hrg = 30;
f = 1;
ref = [17.62040631	0	0	17.62040631
];

[d, h, zone] = test_profile(5);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 12;
path = 3;
htg = 50;
hrg = 10;
f = 2.5;
ref = [25.1016907	0	0	25.1016907
];

[d, h, zone] = test_profile(5);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 13;
path = 3;
htg = 20;
hrg = 20;
f = 0.6;
ref = [20.25488733	0	0	20.25488733

];

[d, h, zone] = test_profile(5);


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

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 14;
path = 3;
htg = 40;
hrg = 50;
f = 0.2;
ref = [10.1899085	0	0	10.1899085
];

[d, h, zone] = test_profile(5);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 15;
path = 3;
htg = 40;
hrg = 50;
f = 0.2;
ref = [10.1899085	0	0	10.1899085
];

[d, h, zone] = test_profile(5);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 16;
path = 3;
htg = 70;
hrg = 5;
f = 0.15;
ref = [16.27593106	6.17409081	10.63765222	20.73949248

];

[d, h, zone] = test_profile(5);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 17;
path = 4;
htg = 30;
hrg = 30;
f = 1;
ref = [7.870985135	0	0	7.870985135
];

[d, h, zone] = test_profile(6);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end


%%
% import test profile

step = 18;
path = 4;
htg = 50;
hrg = 10;
f = 2.5;
ref = [0	0	0	0
];

[d, h, zone] = test_profile(6);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end


%%
% import test profile

step = 19;
path = 4;
htg = 20;
hrg = 20;
f = 0.6;
ref = [15.30876397	0	0	15.30876397
];

[d, h, zone] = test_profile(6);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 20;
path = 4;
htg = 40;
hrg = 50;
f = 0.2;
ref = [6.526730443	0	0	6.526730443
];

[d, h, zone] = test_profile(6);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

%%
% import test profile

step = 21;
path = 4;
htg = 70;
hrg = 5;
f = 0.15;
ref = [0	0	0	0

];

[d, h, zone] = test_profile(6);


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

[hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, htg, hrg, ae, f);

dtot = d(end)-d(1);

%Tx and Rx antenna heights above mean sea level amsl (m)
hts = h(1) + htg;
hrs = h(end) + hrg;

% Compute the path fraction over see

omega = path_fraction(d, zone, 3);

Lbulla = dl_bull(d, h, hts, hrs, ae, f);

% Use the method in 4.2.1 for a second time, with all profile heights hi
% set to zero and modified antenna heights given by

hts1 = hts - hstd;   % eq (38a)
hrs1 = hrs - hsrd;   % eq (38b)
h1 = zeros(size(h));

% where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
% resulting Bullington diffraction loss for this smooth path to Lbulls

Lbulls = dl_bull(d, h1, hts1, hrs1, ae, f);

% Use the method in 4.2.2 to calculate the spherical-Earth diffraction loss
% for the actual path length (dtot) with 

hte = hts1;             % eq (39a)
hre = hrs1;             % eq (39b)
dtot = d(end) - d(1);

Ldsph = dl_se(dtot, hte, hre, ae, f, omega);

% Diffraction loss for the general paht is now given by

Ld(1) = Lbulla + max(Ldsph(1) - Lbulls, 0);  % eq (40)
Ld(2) = Lbulla + max(Ldsph(2) - Lbulls, 0);  % eq (40)%%



comp = [Lbulla(1), Lbulls(1), Ldsph(1), Ld(1)];

error = max(abs(comp-ref));

error = max(abs(comp-ref));

if (error < 0.5)
    fprintf('%d...  Path %d, htg = %.1f, hrg = %.1f f = %.1f passed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    succ = succ + 1;
else
    fprintf(1,'%d...  Path %d, htg = %f, hrg = %f f = %.1f failed\n',step, path, htg, hrg,f)
    fprintf(1,'     (Maximum distance from reference (dB): %f)\n', error);
    fail = fail + 1;
end

return
end



