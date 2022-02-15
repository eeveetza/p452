% MATLAB/Octave script that is used to verify the implementation of
% Recommendation ITU-R P.452-16 (as defined in the file tl_p452.m and the
% functions called therefrom) using a set of test terrain profiles provided by the user.
%
% The script reads all the test profiles from the folder defined by
% the variable <test_profiles>, calculates path profile parameters and
% compares the computed basic transmission loss with the reference ones.

% Author: Ivica Stevanovic (IS), Federal Office of Communications, Switzerland
% Revision History:
% Date            Revision
% 13JUL21         Renaming subfolder "src" into "private" which is automatically in the MATLAB search path
%                                                       (as suggested by K. Konstantinou, Ofcom UK)
% 05JUN2020       Introduced Octave specific code (with M. Rohner, LS telcom)
% 16JUN2016       Initial version (IS)


if (ispc)
    
    clear all;
    close all;
    fclose all;
    
    tol = 1e-4;
    success = 0;
    total = 0;
    
    %% compute the path profile parameters
    s = pwd;
    % if ~exist('p676d11_ga.m','file')
    %     addpath ([s '/src/'])
    % end
    
    % path to the folder containing test profiles
    test_profiles = [s '/validation_examples/'];
    
    
    
    %% begin code
    % Collect all the filenames .csv in the folder pathname that contain the profile data
    filenames = dir(fullfile(test_profiles, '*.xlsx')); % filenames(i).name is the filename
    
    try
        % start excel application
        if(isOctave)
            pkg load windows
        end
        e = actxserver ('Excel.Application');
        %for iname = 1:1
        for iname = 1 : length(filenames)
            %filename1 = 'test_profile_mixed_109km.xlsx'
            filename1 = filenames(iname).name;
            fprintf(1,'***********************************************\n');
            fprintf(1,'Processing file %s ...\n', filename1);
            fprintf(1,'***********************************************\n');
            
            
            failed = false;
            clear p452 z
            
            
            %% open the excel file
            if (~isOctave) %Matlab specific code for reading Excel files
                fid = e.Workbooks.Open([test_profiles filename1]);
                %% read the worksheet path profile
                shid = e.Worksheets.get('Item','Path profile');
                n =  shid.Range('A2').End('xlDown').Row;
                col = shid.Range(['A2:A' num2str(n)]).Value;
                p452.path.d = cell2mat(col).';
                col = shid.Range(['B2:B' num2str(n)]).Value;
                p452.path.h = cell2mat(col).';
                col = shid.Range(['C2:C' num2str(n)]).Value;
                for i = 1:length(col)
                    if strcmp(col{i}, 'A2')
                        z(i) = 2;
                    elseif strcmp(col{i}, 'A1')
                        z(i) = 1;
                    elseif strcmp(col{i}, 'B')
                        z(i) = 3;
                    end
                end
                p452.path.zone = z;
                
                %% read the worksheet Input parameters
                shid = e.Worksheets.get('Item','Input parameters');
                n =  shid.Range('A2').End('xlDown').Row;
                instr = shid.Range(['A2:A' num2str(n)]).Value;
                value = shid.Range(['B2:B' num2str(n)]).Value;
                for i = 1: length(value)
                    p452.(instr{i}) = value{i};
                end
                if p452.polarization == 'v'
                    p452.pol = 2;
                else
                    p452.pol = 1;
                end
                
                %% read reference path profile parameters
                
                shid = e.Worksheets.get('Item','Path profile parameters');
                n =  shid.Range('A2').End('xlDown').Row;
                instr = shid.Range(['A2:A' num2str(n)]).Value;
                value = shid.Range(['B2:B' num2str(n)]).Value;
                for i = 1: length(value)
                    ppref.(instr{i}) = value{i};
                end
                
                if strcmp(ppref.path,'Line of Sight')
                    ppref.pathtype = 1;
                else
                    ppref.pathtype = 2;
                end
                
            else % Octave specific code for reading Excel files
                
                %% open the excel file
                
                fid = e.Workbooks.Open([test_profiles filename1]);
                %% read the worksheet path profile
                shid = e.Worksheets(1);
                
                ur = shid.UsedRange.Value;
                %n =  size(ur)(1)
                n =  size(ur,1);
                
                col = ur(2:n,1);
                p452.path.d = cell2mat(col).';
                col = ur(2:n,2);
                p452.path.h = cell2mat(col).';
                col = ur(2:n,3);
                col(strcmp ("A2", col)) = 2;
                col(strcmp ("A1", col)) = 1;
                col(strcmp ("B", col)) = 3;
                p452.path.zone = cell2mat(col).';
                
                %% read the worksheet Input parameters
                shid = e.Worksheets(2);
                ur = shid.UsedRange.Value;
                %n = size(ur(:,1))(1);
                n = size(ur(:,1),1);
                instr = shid.Range(['A2:A' num2str(n)]).Value;
                value = shid.Range(['B2:B' num2str(n)]).Value;
                for i = 1: length(value)
                    p452.(instr{i}) = value{i};
                end
                if p452.polarization == 'v'
                    p452.pol = 2;
                else
                    p452.pol = 1;
                end
                
                %% read reference path profile parameters
                
                shid = e.Worksheets(3);
                ur = shid.UsedRange.Value;
                n = size(ur(:,1),1);
                instr = shid.Range(['A2:A' num2str(n)]).Value;
                value = shid.Range(['B2:B' num2str(n)]).Value;
                for i = 1: length(value)
                    ppref.(instr{i}) = value{i};
                end
                
                if strcmp(ppref.path,'Line of Sight')
                    ppref.pathtype = 1;
                else
                    ppref.pathtype = 2;
                end
                
            end
            
            %[dc, hc, zonec, htgc, hrgc, Aht, Ahr] = closs_corr(p452.f, p452.path.d, p452.path.h, p452.path.zone, p452.htg, p452.hrg, p452.ha_t, p452.ha_r, p452.dk_t, p452.dk_r);
            
            dc = p452.path.d;
            hc = p452.path.h;
            zonec = p452.path.zone;
            htgc = p452.htg;
            hrgc = p452.hrg;
            
            % Path center latitude
            % (great circle calculation according to
            % P.2001 Annex H maybe more appropriate for longer paths)
            phi_path = (p452.phi_t + p452.phi_r)/2;
            
            % Compute  dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
            zone_r = 12;
            dtm = longest_cont_dist(p452.path.d, p452.path.zone, zone_r);
            
            % Compute  dlm     -   the longest continuous inland section of the great-circle path (km)
            zone_r = 2;
            dlm = longest_cont_dist(p452.path.d, p452.path.zone, zone_r);
            
            % Compute b0
            b0 = beta0(phi_path, dtm, dlm);
            
            [ae, ab] = earth_rad_eff(p452.DN);
            
            [hst, hsr, hstd, hsrd, hte, hre, hm, dlt, dlr, theta_t, theta_r, theta, pathtype] = smooth_earth_heights(dc, hc, htgc, hrgc, ae, p452.f);
            
            dtot = dc(end)-dc(1);
            
            %Tx and Rx antenna heights above mean sea level amsl (m)
            hts = hc(1) + htgc;
            hrs = hc(end) + hrgc;
            
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
            
            %% verify the results out against reference ppref
            if (~isOctave)
                flds = fields(out);
                for i = 1:length(flds)
                    error = abs(out.(flds{i})-ppref.(flds{i}));
                    if error > tol
                        fprintf(1,'Error in %s larger than tolerance %g: %g\n', flds{i}, tol, error);
                        failed = true;
                    end
                end
            end
            
            %% read reference transmission loss parameters
            if (~isOctave) % Matlab specific functions for reading Excel
                shid = e.Worksheets.get('Item','Transmission losses');
                n =  shid.Range('A3').End('xlDown').Row;
                ff = shid.Range(['A3:A' num2str(n)]).Value;
                pp = shid.Range(['B3:B' num2str(n)]).Value;
                Lb_ref = shid.Range(['C3:C' num2str(n)]).Value;
                Lbfsg_ref = shid.Range(['D3:D' num2str(n)]).Value;
                Lb0p_ref = shid.Range(['E3:E' num2str(n)]).Value;
                Lb0b_ref = shid.Range(['F3:F' num2str(n)]).Value;
                Ldsph_ref = shid.Range(['G3:G' num2str(n)]).Value;
                Ld50_ref = shid.Range(['H3:H' num2str(n)]).Value;
                Ldp_ref = shid.Range(['I3:I' num2str(n)]).Value;
                Lbs_ref = shid.Range(['J3:J' num2str(n)]).Value;
                Lba_ref = shid.Range(['K3:K' num2str(n)]).Value;
                
            else % octave specific functions for reading Excel files
                
                shid = e.Worksheets(4);
                ur = shid.UsedRange.Value;
                %n = size(ur(:,1))(1);
                n = size(ur(:,1),1);
                ff = shid.Range(['A3:A' num2str(n)]).Value;
                pp = shid.Range(['B3:B' num2str(n)]).Value;
                Lb_ref = shid.Range(['C3:C' num2str(n)]).Value;
                Lbfsg_ref = shid.Range(['D3:D' num2str(n)]).Value;
                Lb0p_ref = shid.Range(['E3:E' num2str(n)]).Value;
                Lb0b_ref = shid.Range(['F3:F' num2str(n)]).Value;
                Ldsph_ref = shid.Range(['G3:G' num2str(n)]).Value;
                Ld50_ref = shid.Range(['H3:H' num2str(n)]).Value;
                Ldp_ref = shid.Range(['I3:I' num2str(n)]).Value;
                Lbs_ref = shid.Range(['J3:J' num2str(n)]).Value;
                Lba_ref = shid.Range(['K3:K' num2str(n)]).Value;
                
            end
            % compute the transmission losses from MATLAB functions
            
            offset = 0;
            for i = 1:n-2
                [Lbfsg{offset + i}, Lb0p{offset + i}, Lb0b{offset + i}] = pl_los(dtot, ...
                    ff{i}, ...
                    pp{i}, ...
                    b0, ...
                    omega, ...
                    p452.temp, ...
                    p452.press,...
                    dlt, ...
                    dlr);
                
                Lbs{offset + i} = tl_tropo(dtot, ...
                    theta, ...
                    ff{i}, ...
                    pp{i}, ...
                    p452.temp, ...
                    p452.press, ...
                    p452.N0, ...
                    p452.Gt, ...
                    p452.Gr );
                
                
                Lba{offset + i} = tl_anomalous(dtot, ...
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
                    ff{i}, ...
                    pp{i}, ...
                    p452.temp, ...
                    p452.press, ...
                    omega, ...
                    ae, ...
                    b0);
                
                Lbulla{offset + i} = dl_bull(dc, hc, hts, hrs, ae, ff{i});
                
                % Use the method in 4.2.1 for a second time, with all profile heights hi
                % set to zero and modified antenna heights given by
                
                hts1 = hts - hstd;   % eq (38a)
                hrs1 = hrs - hsrd;   % eq (38b)
                h1 = zeros(size(hc));
                
                % where hstd and hsrd are given in 5.1.6.3 of Attachment 2. Set the
                % resulting Bullington diffraction loss for this smooth path to Lbulls
                
                Lbulls{offset + i} = dl_bull(dc, h1, hts1, hrs1, ae, ff{i});
                
                % Use the method in 4.2.2 to radiomaps the spherical-Earth diffraction loss
                % for the actual path length (dtot) with
                
                hte1 = hts1;             % eq (39a)
                hre1 = hrs1;             % eq (39b)
                
                Ldsph{offset + i} = dl_se(dtot, hte1, hre1, ae, ff{i}, omega);
                
                % Diffraction loss for the general paht is now given by
                
                Ld{offset + i}(1) = Lbulla{offset + i} + max(Ldsph{offset + i}(1) - Lbulls{offset + i}, 0);  % eq (40)
                Ld{offset + i}(2) = Lbulla{offset + i} + max(Ldsph{offset + i}(2) - Lbulls{offset + i}, 0);  % eq (40)%%
                
                [ Ldp{offset+i}, Ld50{offset+i} ] = dl_p( dc, hc, hts, hrs, hstd, hsrd, ff{i}, omega, pp{i}, b0, p452.DN );
                
                Lb{offset+i} = tl_p452_pdr(ff{i}, ...
                    pp{i}, ...
                    p452.path.d, ...
                    p452.path.h, ...
                    p452.path.zone, ...
                    p452.path.h, ... % clutter + terrain profile along the path - here we assume clutter = 0 for validation purposes
                    p452.htg, ...
                    p452.hrg, ...
                    p452.phi_t,...
                    p452.phi_r, ...
                    p452.Gt, ...
                    p452.Gr, ...
                    p452.pol, ...
                    p452.dct, ...
                    p452.dcr, ...
                    p452.DN, ...
                    p452.N0, ...
                    p452.press, ...
                    p452.temp);
                
                out1.Lbfsg(i+offset) = Lbfsg{i}-Lbfsg_ref{i};
                out1.Lb0p(i+offset) = Lb0p{i}-Lb0p_ref{i};
                out1.Lb0b(i+offset) = Lb0b{i}-Lb0b_ref{i};
                out1.Ldsph(i+offset) = Ldsph{i}(p452.pol)-Ldsph_ref{i};
                out1.Ld50(i+offset) = Ld50{i}(p452.pol)-Ld50_ref{i};
                out1.Ldp(i+offset) = Ldp{i}(p452.pol)-Ldp_ref{i};
                out1.Lbs(i+offset) = Lbs{i}-Lbs_ref{i};
                out1.Lba(i+offset) = Lba{i}-Lba_ref{i};
                out1.Lb(i+offset) = Lb{i}-Lb_ref{i};
                
                
            end
            
            %% verify error in results out1 against tolarance
            if (~isOctave)
                flds = fields(out1);
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
            end
            
            
            if(isOctave)
                delete (shid);
            end
            fid.Close(false);
            
            if (~failed)
                success = success + 1;
            end
            total = total + 1;
            
        end
        
        %# close Excel
        if(isOctave)
            e.Quit();
            delete(e);
        else
            e.Quit();
            e.delete();
        end
        
        fprintf(1, 'Validation results: %d out of %d tests passed successfully.\n', success, total);
        if (success == total)
            fprintf(1,'The deviation from the reference results is smaller than %g dB.\n', tol);
        end
        
    catch ME
        if(isOctave)
            e.Quit();
            delete(e);
        else
            e.Quit();
            e.delete();
        end
        rethrow(ME)
    end
    
else
    disp('Platform not supported. This script runs only on Windows.')
end
