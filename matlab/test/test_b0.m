% Function being tested: beta_0

disp('Testing function beta0');

% add path to the folder where the functions are defined
s = pwd;
s=s(1:end-5);
if (~exist('longest_cont_dist.m','file'))
    addpath([s '/src'])
end
if (~exist('tl_p452.m','file'))
    addpath(s)
end

% import test profile

[d, h, zone] = test_profile(1);

% Latitude of Tx station
phi_t = 51.2;

% Latitude of Rx station
phi_r = 50.73;

% Path center latitude

phi_path = (phi_t + phi_r)/2;

% Meaning of zone
% 1 = A1 Coastal land: Coastal land and shore areas, i.e. land adjacent to the sea up to an altitude of 100 m 							
% 		relative to mean sea or water level, but limited to a distance of 50 km from the nearest 							
% 		sea area. Where precise 100 m data are not available an approximate value, i.e. 							
% 		300 ft, may be used							
% 2 = A2    Inland:	All land, other than coastal and shore areas defined as “coastal land” above							
% 3 = B		Seas, oceans and other large bodies of water (i.e. covering a circle of at least 100 km 							
% 		in diameter)	

% Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)

zone_r = 1;

dtm = longest_cont_dist(d, zone, zone_r);

% Compute  dlm     -   the longest continuous inland section of the great-circle path (km)

zone_r = 2;

dlm = longest_cont_dist(d, zone, zone_r);


b0 = beta0(phi_path, dtm, dlm);
b0_itur = 3.2827313890624400;

error = -log10(abs(b0-b0_itur));

if (error > 14)
    disp('... passed')
else
    fprintf(1,'...failed\n');
    fprintf(1, 'b0 = %g\n', b0);
    fprintf(1, 'bref = %g\n', b0_itur);
end


