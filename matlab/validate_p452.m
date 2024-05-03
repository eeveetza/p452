% MATLAB/Octave script that is used to verify the implementation of
% Recommendation ITU-R P.452 (as defined in the file tl_p452.m and the
% functions called therefrom) using a set of test terrain profiles provided by the user.
%
% The script reads all the test profiles from the folder defined by
% the variable <test_profiles>, calculates path profile parameters and
% compares the computed basic transmission loss with the reference ones.

% Author: Ivica Stevanovic (IS), Federal Office of Communications, Switzerland
% Revision History:
% Date            Revision
% 15NOV23         Updated to align with ITU-R P.452-18 (IS)
% 16FEB22         Created - transliterated validation examples using .csv instead of .xlsx
%                 Created new validate_p452.m that handles .csv files
%                 These are meant to replace the previous version working
%                 with Excel files. The advantage: it works in both MATLAB
%                 and Octave, Windows and MacOS and it is easier to diff
% 13JUL21         Renaming subfolder "src" into "private" which is automatically in the MATLAB search path
%                                                       (as suggested by K. Konstantinou, Ofcom UK)
% 05JUN2020       Introduced Octave specific code (with M. Rohner, LS telcom)
% 16JUN2016       Initial version (IS)


clear all;
close all;
fclose all;

tol = 1e-6;
success = 0;
total = 0;

flag_createlog = 0;

%% compute the path profile parameters
s = pwd;

% path to the folder containing test profiles
test_profiles = [s '/validation_examples/profiles/'];
test_results  = [s '/validation_examples/results/'];



%% begin code
% Collect all the filenames .csv in the folder pathname that contain the profile data
filenames = dir(fullfile(test_profiles, '*.csv')); % filenames(i).name is the filename


% start excel application
if(isOctave)
    pkg load windows
end

for iname = 1 : length(filenames)
    
    filename1 = filenames(iname).name;
    fprintf(1,'***********************************************\n');
    fprintf(1,'Processing file %s ...\n', filename1);
    fprintf(1,'***********************************************\n');
    
    
    failed = false;
    clear p452 z
    
    
    % read the path profile
    X = readcsv([test_profiles  filename1]);
    
    p452.path.d = str2double( X(:,1) );
    p452.path.h = str2double( X(:,2) );
    p452.path.r = str2double( X(:,3) );
    p452.path.zone = str2double( X(:,5) );
    
    p452.path.g = p452.path.h + p452.path.r;
   
    % Apply the condition in Step 4: Radio profile 
    % gi is the terrain height in metres above sea level for all the points at a distance from transmitter or receiver less than 50 m.
    
    kk = find(p452.path.d < 50/1000);
    if (~isempty(kk))
        p452.path.g(kk) = p452.path.h(kk);
    end
    
    kk = find(p452.path.d(end)-p452.path.d < 50/1000);
    if (~isempty(kk))
        p452.path.g(kk) = p452.path.h(kk);
    end
    
    fname_part = filename1(13:end);
    test_result = ['test_result' fname_part] ;
    if(flag_createlog)
        fidlog = fopen(test_result, 'w');
        fprintf(fidlog,'profile,f (GHz),p (%%),htg (m),hrg (m),phit_e (deg),phit_n (deg),phir_e (deg),phir_n (deg),Gt (dBi),Gr (dBi),pol (1-h/2-v),dct (km),dcr (km),press (hPa),temp (deg C),ae,dtot,hts,hrs,theta_t,theta_r,theta,hm,hte,hre,hstd,hsrd,dlt,dlr,path,dtm,dlm,b0,omega,DN,N0,Lb,Lbfsg,Lb0p,Lb0b,Ldsph,Ld50,Ldp,Lbs,Lba\n');
    end
    % read the input arguments and reference values
    
    Y = readcsv([test_results  test_result]);
    
    
    [nrows, ncols] = size(Y);
    ff = zeros(nrows,1);
    for i = 1:nrows
        ff(i)    = str2double(Y(i,2));
        pp(i)    = str2double(Y(i,3));
    end
    
    p452.htg      = str2double(Y(1,4));
    p452.hrg      = str2double(Y(1,5));
    p452.phit_e   = str2double(Y(1,6));
    p452.phit_n   = str2double(Y(1,7));
    p452.phir_e   = str2double(Y(1,8));
    p452.phir_n   = str2double(Y(1,9));

    p452.Gt      = str2double(Y(1,10));
    p452.Gr      = str2double(Y(1,11));
    p452.pol     = str2double(Y(1,12));
    p452.dct     = str2double(Y(1,13));
    p452.dcr     = str2double(Y(1,14));


    p452.press   = str2double(Y(1,15));
    p452.temp    = str2double(Y(1,16));

    
    ppref.ae      =  str2double(Y(1,17));
    ppref.dtot    =  str2double(Y(1,18));
    ppref.hts     =  str2double(Y(1,19));
    ppref.hrs     =  str2double(Y(1,20));
    ppref.theta_t =  str2double(Y(1,21));
    ppref.theta_r =  str2double(Y(1,22));
    ppref.theta   =  str2double(Y(1,23));
    ppref.hm      =  str2double(Y(1,24));
    ppref.hte     =  str2double(Y(1,25));
    ppref.hre     =  str2double(Y(1,26));
    ppref.hstd    =  str2double(Y(1,27));
    ppref.hsrd    =  str2double(Y(1,28));
    ppref.dlt     =  str2double(Y(1,29));
    ppref.dlr     =  str2double(Y(1,30));
    ppref.path    =  (Y(1,31));
    ppref.dtm     =  str2double(Y(1,32));
    ppref.dlm     =  str2double(Y(1,33));
    ppref.b0      =  str2double(Y(1,34));
    ppref.omega   =  str2double(Y(1,35));
    ppref.DN      =  str2double(Y(1,36));
    ppref.N0      =  str2double(Y(1,37));
    
    if strcmp(ppref.path,'Line of Sight')
        ppref.pathtype = 1;
    else
        ppref.pathtype = 2;
    end
    
    d = p452.path.d;
    h = p452.path.h;
    g = p452.path.g;
    zone = p452.path.zone;
    
    % Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
    zone_r = 12;
    dtm = longest_cont_dist(p452.path.d, p452.path.zone, zone_r);
    
    % Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
    zone_r = 2;
    dlm = longest_cont_dist(p452.path.d, p452.path.zone, zone_r);

    % Calculate the longitude and latitude of the mid-point of the path, phim_e,
    % and phim_n for dpnt = 0.5dt
    Re = 6371;
    dpnt = 0.5*(d(end)-d(1));
    [phim_e, phim_n, bt2r, dgc] = great_circle_path(p452.phir_e, p452.phit_e, p452.phir_n, p452.phit_n, Re, dpnt);


    % Find radio-refractivity lapse rate dN
    % using the digital maps at phim_e (lon), phim_n (lat) - as a bilinear interpolation

    DN = get_interp2('DN50',phim_e,phim_n);
    N0 = get_interp2('N050',phim_e,phim_n);


    % Compute b0
    b0 = beta0(phim_n, dtm, dlm);
    
    [ae, ab] = earth_rad_eff(DN);
    
    [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(d, h, p452.htg, p452.hrg, ae, ff(1));
    
    dtot = d(end)-d(1);
    
    %Tx and Rx antenna heights above mean sea level amsl (m)
    hts = h(1) + p452.htg;
    hrs = h(end) + p452.hrg;
    
    % Compute the path fraction over see
    
    omega = path_fraction(p452.path.d, p452.path.zone, 3);
    
    out.ae = ae;
    out.dtot = dtot;
    out.hts = hts;
    out.hrs = hrs;
    out.theta_t = theta_t;
    out.theta_r = theta_r;
    out.theta = theta;
    out.hm = hm;
    out.hte = hte;
    out.hre = hre;
    out.hstd = hstd;
    out.hsrd = hsrd;
    out.dlt = dlt;
    out.dlr = dlr;
    out.pathtype = pathtype;
    out.dtm = dtm;
    out.dlm = dlm;
    out.b0 = b0;
    out.omega = omega;
    out.DN = DN;
    out.N0 = N0;

    %% verify the results struct `out` against the reference struct `ppref`
    
    flds = fieldnames(out);
    for i = 1:length(flds)
        error = abs(out.(flds{i})-ppref.(flds{i}));
        if error > tol
            fprintf(1,'Error in %s larger than tolerance %g: %g\n', flds{i}, tol, error);
            failed = true;
        end
    end
    
    % extract reference transmission losses
    
    Lb_ref    = str2double(Y(:,38));
    Lbfsg_ref = str2double(Y(:,39));
    Lb0p_ref  = str2double(Y(:,40));
    Lb0b_ref  = str2double(Y(:,41));
    Ldsph_ref = str2double(Y(:,42));
    Ld50_ref  = str2double(Y(:,43));
    Ldp_ref   = str2double(Y(:,44));
    Lbs_ref   = str2double(Y(:,45));
    Lba_ref   = str2double(Y(:,46));
    
    
    % compute the transmission losses using MATLAB functions
    Lbfsg = zeros(nrows,1);
    Lb0p = zeros(nrows, 1);
    Lb0b = zeros(nrows, 1);
    Lbs = zeros(nrows, 1);
    Lba = zeros(nrows, 1);
    Lbulla = cell(nrows, 1);
    Lbulls = cell(nrows, 1);
    Ldsph = cell(nrows, 1);
    Ld = cell(nrows, 1);
    Ldp = cell(nrows, 1);
    Ld50 = cell(nrows, 1);
    Lb = zeros(nrows,1);
    
    offset = 0;
    d3D = sqrt(dtot*dtot + ((hts-hrs)/1000.0).^2);
    
    for i = 1:nrows
        
        [Lbfsg(offset + i), Lb0p(offset + i), Lb0b(offset + i)] = pl_los(d3D, ...
            ff(i), ...
            pp(i), ...
            b0, ...
            omega, ...
            p452.temp, ...
            p452.press,...
            dlt, ...
            dlr);

%         % The path length expressed as the angle subtended by d km at the center of
%         % a sphere of effective Earth radius ITU-R P.2001-4 (3.5.4)
% 
%         theta_e = dtot/ae; % radians
% 
%         [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(p452.path.d, p452.path.h, p452.htg, p452.hrg, ae, ff(i));
% 
%         % Calculate the horizon elevation angles limited such that they are positive
% 
%         theta_tpos = max(theta_t, 0);                   % Eq (3.7.11a) ITU-R P.2001-4
%         theta_rpos = max(theta_r, 0);                   % Eq (3.7.11b) ITU-R P.2001-4
% 
%         [dt_cv, phi_cve, phi_cvn] = tropospheric_path(dtot, hts, hrs, theta_e, theta_tpos, theta_rpos, p452.phir_e, p452.phit_e, p452.phir_n, p452.phit_n, Re);
% 
%         % height of the Earth's surface above sea level where the common volume is located
% 
%         Hs = surface_altitude_cv(h, d, dt_cv)/1000.0; % in km
% 
%         [Lbs(offset + i), theta_s] = tl_troposcatter(ff(i), dtot, hts, hrs, ae, theta_e, theta_t, theta_r, phi_cvn, phi_cve, p452.Gt, p452.Gr, pp(i), Hs);
% 
%         %% To avoid under-estimating troposcatter for short paths, limit Lbs (E.17)
%         Lbs(offset + i) = max(Lbs(offset + i), Lbfsg(offset + i));

        
        % Calculate the basic transmission loss due to troposcatter not exceeded
        % for any time percantage p 
        
        Lbs(offset + i) = tl_tropo(dtot, theta, ff(i), pp(i), p452.temp, p452.press, N0, p452.Gt, p452.Gr );
        
        Lba(offset + i) = tl_anomalous(dtot, ...
            dlt, ...
            dlr, ...
            p452.dct, ...
            p452.dcr, ...
            dlm, ...
            hts, ...
            hrs, ...
            hte, ...
            hre, ...
            hm, ...
            theta_t, ...
            theta_r, ...
            ff(i), ...
            pp(i), ...
            p452.temp, ...
            p452.press, ...
            omega, ...
            ae, ...
            b0);
        
        Lbulla{offset + i} = dl_bull(d, g, hts, hrs, ae, ff(i));
        
        % Use the method in 4.2.1 for a second time, with all profile heights hi
        % set to zero and modified antenna heights given by
        
        hts1 = hts - hstd;   % eq (38a)
        hrs1 = hrs - hsrd;   % eq (38b)
        h1 = zeros(size(h));
        
        % where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
        % resulting Bullington diffraction loss for this smooth path to Lbulls
        
        Lbulls{offset + i} = dl_bull(d, h1, hts1, hrs1, ae, ff(i));
        
        % Use the method in 4.2.2 to radiomaps the spherical-Earth diffraction loss
        % for the actual path length (dtot) with
        
        hte1 = hts1;             % eq (39a)
        hre1 = hrs1;             % eq (39b)
        
        Ldsph{offset + i} = dl_se(dtot, hte1, hre1, ae, ff(i), omega);
        
        % Diffraction loss for the general paht is now given by
        
        Ld{offset + i}(1) = Lbulla{offset + i} + max(Ldsph{offset + i}(1) - Lbulls{offset + i}, 0);  % eq (40)
        Ld{offset + i}(2) = Lbulla{offset + i} + max(Ldsph{offset + i}(2) - Lbulls{offset + i}, 0);  % eq (40)%%
        
        [ Ldp{offset+i}, Ld50{offset+i} ] = dl_p( d, g, hts, hrs, hstd, hsrd, ff(i), omega, pp(i), b0, DN );
        
        Lb(offset+i) = tl_p452(ff(i), ...
            pp(i), ...
            p452.path.d, ...
            p452.path.h, ...
            p452.path.g,...
            p452.path.zone, ...
            p452.htg, ...
            p452.hrg, ...
            p452.phit_e,...
            p452.phit_n,...
            p452.phir_e,...
            p452.phir_n,...
            p452.Gt, ...
            p452.Gr, ...
            p452.pol, ...
            p452.dct, ...
            p452.dcr, ...
            p452.press, ...
            p452.temp);
        
        out1.Lbfsg(i+offset) = Lbfsg(i)-Lbfsg_ref(i);
        out1.Lb0p(i+offset) = Lb0p(i)-Lb0p_ref(i);
        out1.Lb0b(i+offset) = Lb0b(i)-Lb0b_ref(i);
        out1.Ldsph(i+offset) = Ldsph{i}(p452.pol)-Ldsph_ref(i);
        out1.Ld50(i+offset) = Ld50{i}(p452.pol)-Ld50_ref(i);
        out1.Ldp(i+offset) = Ldp{i}(p452.pol)-Ldp_ref(i);
        out1.Lbs(i+offset) = Lbs(i)-Lbs_ref(i);
        out1.Lba(i+offset) = Lba(i)-Lba_ref(i);
        out1.Lb(i+offset) = Lb(i)-Lb_ref(i);

        if pathtype == 1
            pathtypestr = 'Line of Sight';
        else
            pathtypestr = 'Trans-Horizon';
        end
        
        if(flag_createlog==1)
         fprintf(fidlog, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%s,%f,%f,%f,%f,%f,%f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f \n', ...
             Y{i,1}, Y{i,2},Y{i,3},Y{i,4},Y{i,5},Y{i,6},Y{i,7},Y{i,8},Y{i,9},...
             Y{i,10},Y{i,11},Y{i,12},Y{i,13},Y{i,14},Y{i,15},Y{i,16},...
                ae, ...
                dtot, ...
                hts, ...
                hrs, ...
                theta_t, ...
                theta_r, ...
                theta, ...
                hm, ...
                hte, ...
                hre, ...
                hstd, ...
                hsrd, ...
                dlt, ...
                dlr, ...
                pathtypestr, ...
                dtm, ...
                dlm, ...
                b0, ...
                omega, ...
                DN, ...
                N0, ...
                Lb(i), ...
                Lbfsg(i), ...
                Lb0p(i), ...
                Lb0b(i), ...
                Ldsph{i}(p452.pol), ...
                Ld50{i}(p452.pol), ...
                Ldp{i}(p452.pol), ...
                Lbs(i), ...
                Lba(i) ...
             );
        
        end
    end

    %% verify error in the results out1 against tolarance
    
    flds = fieldnames(out1);
    for i = 1:length(flds)
        maxi = max(abs(out1.(flds{i})));
        kk = find(maxi > tol);
        if ~isempty(kk)
            for kki = 1:length(kk)
                fprintf(1,'Maximum deviation found in %s larger than tolerance %g:  %g\n', flds{i}, tol, maxi(kk(kki)));
                failed = true;
            end
        end
    end
    
    if(flag_createlog)
        fclose(fidlog);
    end
    
    
    if (~failed)
        success = success + 1;
    end
    total = total + 1;
    
end


fprintf(1, 'Validation results: %d out of %d tests passed successfully.\n', success, total);
if (success == total)
    fprintf(1,'The deviation from the reference results is smaller than %g dB.\n', tol);
end




