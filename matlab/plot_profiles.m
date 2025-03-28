
close all
clear all
clc

folder = "C:/Users/U80824876/Meetings/ITUR/SG3/UK_50ProfileData/UK/";
test_folder = "London_5850";
test_file = "5850_523775_182825.csv"; % produces some 15 dB difference
%test_file = "5850_526275_184775.csv"; % gives the same value
%test_file = "3602_549375_357025.csv"; % gives some 10 dB difference, but pdr closer to measurements
fileformat = 'Fryderyk_csv';
%test_folder = "Boston_5850";
%test_file = "5850_538025_342425.csv"; % gives the same value
%test_file = "5850_505375_344875.csv"; % gives 13 db difference

% test_folder = "Merthyr_5850";
% test_file = "5850_281675_201875.csv"; % gives the same value
% test_file = "5850_291575_206875.csv"; % gives 10.7 dB difference
% 
% test_folder = "Nottingham_5850";
% test_file = "5850_437475_347775.csv"; % gives the same value
% test_file = "5850_439425_339375.csv"; % 10.7 dB difference


 test_file = "5850_523575_181175.csv"; % the same value
 test_file = "5850_521675_180875.csv"; % 12 dB difference



sg3db=read_sg3_measurements(strjoin([folder test_folder '/' test_file],''),fileformat);



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

count = 1;
imin = 10;

for i = imin:length(sg3db.x)
    dd = sg3db.x(1:i);
    hh = sg3db.h_gamsl(1:i);
    RR = R(1:i);
    cc = sg3db.coveragecode(1:i);
    zz = climzone(1:i);

    Re = 6371;
    [Phipnte, Phipntn, Bt2r, dgc] = great_circle_path(Phire, Phite, Phirn, Phitn, Re, dd(end));

    [Lb1(i-imin+1), Lbs1(i-imin+1)] = tl_p452_pdr(f/1e3, ...
                    sg3db.TimePercent, ...
                    dd, ... 
                    hh, ...
                    hh + RR, ...
                    zz, ...
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
                    false);

   [Lb2(i-imin+1), Lbs2(i-imin+1)] = tl_p452_pdr(f/1e3, ...
                    sg3db.TimePercent, ...
                    dd, ... 
                    hh, ...
                    hh + RR, ...
                    zz, ...
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
                    true);

end

% print profile

titlestr = replace([test_folder test_file], "_", "\_");

figure
plot(sg3db.x,sg3db.h_gamsl+R,'LineWidth',1,'Color','g', 'LineStyle', '-')
hold on
plot(sg3db.x,sg3db.h_gamsl,'LineWidth',2,'Color','k', 'LineStyle', '-')

xlabel('Distance (km)')
ylabel('Height (m)')
title(titlestr)
grid on

legend( 'Clutter profile', 'Terrain profile','Location','northwest')

        figurename = strjoin([test_folder '_' test_file '_1_profile.png'],'');
        saveas(gcf,figurename)
figure
plot(sg3db.x(imin:end),Lbs1,'LineWidth',1,'Color',[1 0.6 0.6]);
hold on
plot(sg3db.x(imin:end),Lb1,'LineWidth',1,'Color','r', 'LineStyle', '-')
plot(sg3db.x(imin:end),Lbs2,'LineWidth',1,'Color',[0.6 0.6 1]);
plot(sg3db.x(imin:end),Lb2,'LineWidth',1,'Color','b', 'LineStyle', '--')
plot(sg3db.x(end), PL, '*')


xlabel('Distance (km)')
ylabel('Basic transmission loss (dB)')
title(titlestr)
grid on

legend( 'Lbs P.452-18','Lb P.452-18','Lbs PDR P.452-18','Lb PDR P.452-18','Location','northwest')
        figurename = strjoin([test_folder '_' test_file '_2_loss.png'],'');
        saveas(gcf,figurename)