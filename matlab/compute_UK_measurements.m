% This script computes the basic transmission loss according to ITU-R
% P.1812-6 (using tl_p1812.m) for the UK sub-6 GHz measurement data

clear all
close all
clc

warning off

% presenting results with rounding off to 8 digit
rd = 8;

% discard the results with deviation larger than discard
discard = 1000;

%pdr_type = 1;  label = 'PDR P.452^{0\theta + 4}';    %theta^2 + 4
%pdr_type = 2; label = 'PDR P.452^{4\theta + 4}';   %theta^2 + 4theta + 4
pdr_type = 3; label = 'PDR P.452^{7\theta + 4}';   %theta^2 + 7theta + 4

locations = ["Boston", "London", "Merthyr", "Nottingham", "ScarHill", "Southampton", "Stevenage"];
%locations = ["Merthyr", "Nottingham", "ScarHill", "Southampton", "Stevenage"];
frequencies = ["449", "915", "1802", "2695", "3602", "5850"];

folder = "C:/Users/U80824876/Meetings/ITUR/SG3/UK_50ProfileData/UK/";

for loc = 1:length(locations)
    for ff = 1:length(frequencies)
        if (loc == 5 && ff <=2 ) % There is no ScarHill_449 nor ScarHIll_915
            continue
        end

        
        test_profiles_string = [folder locations{loc} "_" frequencies{ff} "/"];
        test_profiles = strjoin(test_profiles_string, '');
        fprintf(1,'Processing  %s\n', test_profiles);


        % Collect all the filenames .csv in the folder pathname that contain the profile data
        filenames = dir(fullfile(test_profiles, '*.csv')); % filenames(i).name is the filename
        N = length(filenames);
        
        count = 1;
        for i = 1:N

            filename1 = filenames(i).name;
            %     fprintf(1,'***********************************************\n');
            %     fprintf(1,'Processing file %s%s ...\n', test_profiles, filename1);
            %     fprintf(1,'***********************************************\n');
            fileformat = 'Fryderyk_csv';
            sg3db=read_sg3_measurements(strjoin([test_profiles filename1],''),fileformat);
            if (length(sg3db.x) < 4) 
                continue
            end

            for pp = 1:2 % once with and once without PDR

                if (pp == 1)
                    pdr = false;
                    %fprintf(1,'First round, applying P.1812-6\n');

                else
                    pdr = true;
                    % fprintf(1,'Second round, applying P.1812-6 PDR\n');
                end

                sg3db.debug = 0;

                % update the data structure with the Tx Power (kW)
                for kindex=1:sg3db.Ndata
                    PERP= sg3db.ERPMaxTotal(kindex);
                    HRED= sg3db.HRPred(kindex);
                    PkW=10^(PERP/10)*1e-3; %kW

                    if(isnan(PkW))
                        % use complementary information from Basic Transmission Loss and
                        % received measured strength to compute the transmitter power +
                        % gain
                        E=sg3db.MeasuredFieldStrength(kindex);
                        PL=sg3db.BasicTransmissionLoss(kindex);
                        f=sg3db.frequency(kindex);
                        PdBkW=-137.2217+E-20*log10(f)+PL;
                        PkW=10^(PdBkW/10);
                    end

                    sg3db.TransmittedPower(kindex)=PkW;
                end



                % transform coverage code to the representative clutter heights
                %sg3db.h_ground_cover = coverage_code2clutter_height(sg3db.coveragecode);
                R = sg3db.h_ground_cover;
                % Check if R entries are empty or NaN. If they are assign
                % the representative clutter heights according to p1058code
                if(isnan(R(1)))
                    R = p1058code_to_p1812ch(sg3db.coveragecode);
                end

                % zone    -   Zone type: Coastal land (3), Inland (4) or Sea (1)
                climzone_p1812 = sg3db.radio_met_code;
                % the climatic zone codes of P1812 are different to the ones of
                % P.452
                climzone = cz_p1812_to_p452(climzone_p1812);
                sg3db.ClutterCode = p1058code_to_p1812code(sg3db.coveragecode, R);
         

                Phire = sg3db.RxLon;
                Phirn = sg3db.RxLat;
                Phite = sg3db.TxLon;
                Phitn = sg3db.TxLat;
                Hrg = sg3db.hRx;
                Htg = sg3db.hTx;
                Grx = 0;
                Gtx = 0;

                FlagVP = sg3db.polHVC;  % vertical polarization
                dct = 500;
                if sg3db.coveragecode(1) == 1
                    dct = 0;
                end
                dcr = 500;
                if sg3db.coveragecode(end) == 1
                    dcr = 0;
                end
                press = 1013;
                temp = 20;

                Lb = tl_p452_pdr(f/1e3, ...
                    sg3db.TimePercent, ...
                    sg3db.x, ... 
                    sg3db.h_gamsl, ...
                    sg3db.h_gamsl + R, ...
                    climzone, ...
                    Htg, ...
                    Hrg, ...
                    Phite, ...
                    Phitn, ...
                    Phire, ...
                    Phirn, ...
                    Gtx, ...
                    Grx, ...
                    FlagVP, ...
                    dct, ...
                    dcr, ...
                    press, ...
                    temp, ...
                    pdr, ...
                    pdr_type);


                %         fprintf(1,'Measured basic tl: %g dB\n', sg3db.BasicTransmissionLoss);
                %         fprintf(1,'Simulated basic tl: %g dB\n', Lb);
                %         fprintf(1,'Error: %g dB\n', Lb - sg3db.BasicTransmissionLoss);
                %         fprintf(1,'\n');

                result{pp}{count}.statno = num2str(filename1(1:end-4));
                result{pp}{count}.d = sg3db.x(end);
                result{pp}{count}.t = sg3db.TimePercent;
                result{pp}{count}.f = sg3db.frequency;
                result{pp}{count}.plm = round(sg3db.BasicTransmissionLoss,rd);
                result{pp}{count}.pls = round(Lb,rd);
                result{pp}{count}.pls_total = round(Lb,rd);
                result{pp}{count}.delta = round(Lb - sg3db.BasicTransmissionLoss,rd);
                result{pp}{count}.deltaf = Lb - sg3db.BasicTransmissionLoss;



            end

            if (i>=floor(N/100) && i<floor(N/100)+1)
                fprintf(1, '1%% ... ');
            elseif (i>=floor(N/10) && i<floor(N/10)+1)
                fprintf(1, '10%% ... ');
            elseif (i>= floor(0.2*N) && i<floor(0.2*N)+1)
                fprintf(1,'20%% ... ');
            elseif (i>= floor(0.3*N) && i<floor(0.3*N)+1)
                fprintf(1,'30%% ... ');
            elseif (i>= floor(0.4*N) && i<floor(0.4*N)+1)
                fprintf(1,'40%% ... ');
            elseif(i>= floor(N/2) && i<floor(N/2)+1)
                fprintf(1,'50%% ... ');
            elseif (i>= floor(0.6*N) && i<floor(0.6*N)+1)
                fprintf(1,'60%% ... ');
            elseif(i>=floor(0.70*N) && i<floor(0.70*N)+1)
                fprintf(1,'70%% ... ');
            elseif (i>= floor(0.8*N) && i<floor(0.8*N)+1)
                fprintf(1,'80%% ... ');
            elseif(i>=floor(0.9*N) && i<floor(0.9*N)+1)
                fprintf(1,'90%% ... ');
            elseif(i>=floor(0.99*N) && i<floor(0.99*N)+1)
                fprintf(1,'99%% ... ');
            elseif(i==N)
                fprintf(1,'100%% \n');

            end

            count = count + 1;

        end

        dummy = split(test_profiles,'/');
        titlestr1 = dummy{end-1};
        titlestr = replace(titlestr1, "_", "\_");

        filename_out = ['Results_UK_P452_PE_' titlestr1 '_' num2str(pdr_type) '.xls'];

        %fprintf(1,'%10s  %10s  %10s %10s %20s  %20s  %20s  %20s  %20s\n','Stat. no.', 'd (km)', 't (%)', 'f (GHz)', 'Measured PL (dB)', 'P.1812 (dB)', 'PDR (dB)', 'PE P.1812 (dB)', 'PE PDR (dB)');

        A = {'Stat. no.', 't (%)', 'Meas. BTL (dB)', 'P.452 (dB)', 'PDR (dB)', 'PE P.452 (dB)', 'PE PDR (dB)'};

        for kk = 1:length(result{1})
            %fprintf(1,'%10d  %10g  %10g  %10g  %20g  %20g  %20g  %20g  %20g\n', result{1}{kk}.statno, result{1}{kk}.d, result{1}{kk}.t, result{1}{kk}.f, result{1}{kk}.plm, result{1}{kk}.pls_total, result{2}{kk}.pls_total, result{1}{kk}.delta, result{2}{kk}.delta);

            row = {result{1}{kk}.statno, result{1}{kk}.t, result{1}{kk}.plm, result{1}{kk}.pls_total, result{2}{kk}.pls_total, result{1}{kk}.delta, result{2}{kk}.delta};
            if (abs(result{1}{kk}.delta) < discard && abs(result{2}{kk}.delta) < discard)
                A = [A; row];
            end

        end



        if exist(filename_out,'file')

            % if the file already exist, delete it
            [status, result1] = system(['del ' filename_out]);
            fprintf(1,'Rewriting the existing file: %s\n', filename_out);

        end

        fprintf(1,'Writing the results in file: %s\n', filename_out);

        clear B

        B = A(1,:);
        B = [B; A(2:end, :)];

        writecell(B,  filename_out, 'Sheet', 'Sheet1');

        close all
        clear result A B


    end
end
