% This script computes the statistics over all measurements sites
% compute_UK_measurements needs to be run before this script
clear all
close all
clc

warning off

% presenting results with rounding off to 1 digit
rd = 1;

% discard the results with deviation larger than discard
discard = 1000;

pdr_type = 1;  label = 'PDR P.452^{0\theta + 4}';    %theta^2 + 4
%pdr_type = 2; label = 'PDR P.452^{4\theta + 4}';   %theta^2 + 4theta + 4
%pdr_type = 3; label = 'PDR P.452^{7\theta + 4}';   %theta^2 + 7theta + 4

locations = ["Boston", "London", "Merthyr", "Nottingham", "ScarHill", "Southampton", "Stevenage"];
%locations = ["Merthyr", "Nottingham", "ScarHill", "Southampton", "Stevenage"];
frequencies = ["449", "915", "1802", "2695", "3602", "5850"];

folder = "./";

    file_total = ["Results_UK_P452_all_locations_" num2str(pdr_type), ".xls"];
    filename_total = strjoin(file_total, ''); 

    if exist(filename_total,'file')

        % if the file already exist, delete it
        [status, result1] = system(strjoin(['del ' filename_total]),'');
        fprintf(1,'Rewriting the existing file: %s\n', filename_total);

    end

for ff = 1:length(frequencies)

    file_total_freq = ["Results_UK_P452_all_locations_" frequencies{ff} ".xls" ];
    filename_total_freq = strjoin(file_total_freq, ''); 
    titlestr = num2str(frequencies{ff});

    m = [];
    p1 = [];
    p2 = [];


    for loc = 1:length(locations)

        if (loc == 5 && ff <=2 ) % There is no ScarHill_449 nor ScarHIll_915
            continue
        end


        test_profiles_string = [folder locations{loc} "_" frequencies{ff} "/"];
        test_profiles = strjoin(test_profiles_string, '');
        fprintf(1,'Processing  %s\n', test_profiles);

        dummy = split(test_profiles,'/');
        titlestr1 = dummy{end-1};
        titlestr = replace(titlestr1, "_", "\_");

        filename_out = ['Results_UK_P452_PE_' titlestr1 '_' num2str(pdr_type) '.xls'];

        if ~exist(filename_out,'file')
            errormsg = [filename_out " does not exist. Stopping."];
            error(errormsg);
        end


        m_cell  = readcell(filename_out, "NumHeaderLines", 1, "Range", "C:C");
        p1_cell = readcell(filename_out, "NumHeaderLines", 1, "Range", "D:D");
        p2_cell = readcell(filename_out, "NumHeaderLines", 1, "Range", "E:E");

        m_vec  = cell2mat(m_cell(1:end));
        p1_vec = cell2mat(p1_cell(1:end));
        p2_vec = cell2mat(p2_cell(1:end));


        m = [m; m_vec];
        p1 = [p1; p1_vec];
        p2 = [p2; p2_vec];


    end

    delta1 = p1 - m;
    delta2 = p2 - m;

    [max1, min1, muabs1, mu1, median1, sigma1, rmsd1, rho1] = predictionErrorStatistics(m,p1);
    [max2, min2, muabs2, mu2, median2, sigma2, rmsd2, rho2] = predictionErrorStatistics(m,p2);
    
    row = {frequencies{ff}, round(median1,rd), round(mu1,rd), round(sigma1,rd), round(min1,rd), round(max1,rd), round(rho1,rd), length(p1)};
    B = [row];

    writecell(B,  filename_total, 'Sheet', 'In-Force', 'WriteMode', 'append');
    
    row = {frequencies{ff}, round(median2,rd), round(mu2,rd), round(sigma2,rd), round(min2,rd), round(max2,rd), round(rho2,rd), length(p2)};
    B = [row];
    writecell(B,  filename_total, 'Sheet', 'PDR', 'WriteMode', 'append');

    row = {frequencies{ff}, round(mu1,rd), round(mu2,rd), round(sigma1,rd), round(sigma2,rd)};
    B = [row];
    writecell(B,  filename_total, 'Sheet', 'In-force - PDR', 'WriteMode', 'append');
    

    delta_pdr = p2-m;
    delta_inf = p1-m;

    pd_pdr = fitdist(delta_pdr, 'normal');
    pd_inf = fitdist(delta_inf, 'normal');

    figure

    histogram(delta_inf, 'Normalization', 'pdf', 'FaceColor', [1, 1, 0]);
    hold on
    x = [-50:1:80];
    line(x,pdf(pd_inf,x),'LineStyle','-','Color','r');
    lstr1 = ['N(' num2str(pd_inf.mu) ', ' num2str(pd_inf.sigma) ')' ];

    grid on

    histogram(delta_pdr, 'Normalization', 'pdf', 'FaceColor', [1, 0.7, 0]);
    hold on
    x = [-50:1:80];
    line(x,pdf(pd_pdr,x),'LineStyle','-.','Color','b');
    lstr2 = ['N(' num2str(pd_pdr.mu) ', ' num2str(pd_pdr.sigma) ')' ];
    legend('P.452-18', lstr1, label, lstr2)
    xlabel('PE (dB)')
    ylabel('n.u.')
    grid on
    titlestr = num2str(frequencies{ff});
    title(titlestr);
    dummyfigurename = [extractBefore(filename_total_freq, ".xls") '_01_pdf.png'];
    figurename = strjoin(dummyfigurename, ''); 
    saveas(gcf,figurename)


    L_pdr = p2;
    L_inf = p1;
    L_m =   m;

    figure
    plot(L_m, L_inf, 'ro');
    hold on
    plot(L_m, L_pdr,'b.')
    plot(L_m, L_m, 'k')
    xlabel('Lb measured (dB)')
    ylabel('Lb simulated (dB)')
    legend( 'P.452-18',label,'measurements', 'Location','southeast')

    title(titlestr)
    grid on
    
    dummyfigurename = [extractBefore(filename_total_freq, ".xls") '_02_scatt.png'];
    figurename = strjoin(dummyfigurename, ''); 
    saveas(gcf,figurename)

        

    %% prediction error statistics for path loss
    ltmin=100;
    ltmax=180;
    Nlt=1000;
    lt=linspace(ltmin,ltmax,Nlt);
    thr = THRCalc(m, p1,lt);

    kk=find(thr<1);
    if(isempty(kk))
        error('No valid data to process. Program ending...');
    end

    k1=max(kk(1)-1, 1);
    k2=min(kk(end)+1, length(thr));

    % AHRE is the average total hit rate error, large value of AHRE indicates a
    % good fit between the predicted and estimated values
    AHRE=0;
    for k=1:length(kk)
        AHRE=AHRE+thr(kk(k));
    end
    AHRE=AHRE/length(kk)*100;
    kindex=max(1,floor(kk(1)*0.8)):min(floor(k2*1.2),k2);
    figure
    plot(lt(kindex),thr(kindex)*100,'LineWidth',1,'Color','r', 'LineStyle', '--')
    xlabel('Transmission Loss Threshold [dB]')
    ylabel('Total Hit Rate [%]')
    title(titlestr)
    grid on


    hold on


    lt=linspace(ltmin,ltmax,Nlt);
    thr = THRCalc(m, p2,lt);

    kk=find(thr<1);
    if(isempty(kk))
        error('No valid data to process. Program ending...');
    end

    k1=max(kk(1)-1, 1);
    k2=min(kk(end)+1, length(thr));

    % AHRE is the average total hit rate error, large value of AHRE indicates a
    % good fit between the predicted and estimated values
    AHRE=0;
    for k=1:length(kk)
        AHRE=AHRE+thr(kk(k));
    end
    AHRE=AHRE/length(kk)*100;
    kindex=max(1,floor(kk(1)*0.8)):min(floor(k2*1.2),k2);
    hold on
    plot(lt(kindex),thr(kindex)*100,'LineWidth',1,'Color','b')


    legend( 'P.452-18',label,'Location','southwest')
    dummyfigurename = [extractBefore(filename_total_freq, ".xls") '_03_ahre.png'];
    figurename = strjoin(dummyfigurename, ''); 
    
    saveas(gcf,figurename)



end
close all
