% Function being tested: earth_rad_eff

disp('Testing function earth_rad_eff');

% add path to the folder where the functions are defined
s = pwd;
s=s(1:end-5);
if (~exist('longest_cont_dist.m','file'))
    addpath([s '/src'])
end
if (~exist('tl_p452.m','file'))
    addpath(s)
end

DN = 53;

[ae, ab] = earth_rad_eff(DN);

ae_itur = 9617.75961538462000000000;


error = -log10(abs(ae-ae_itur));

if (error > 11)
    disp('... passed')
else
    fprintf(1,'...failed\n');
    fprintf(1, 'ae = %g\n', ae);
    fprintf(1, 'ae_ref = %g\n', ae_itur);
end


