%% MATLAB script used to analyze the behavior of PDR P.452-clutter
% For selected corner cases
% Scenario 1: Flat terrain
% flat = true
% dist = 10, 30, 60 m (free-space, spherical-earth diffraction, tropospheric and other effects)
% Scenario 2: Ideal Triangular "Mountain"
% flat = false
% dist = 10, 30, 60
% The script computes and plots
% Transmission loss as a function of distance (w/o and w/ clutter)
% Excess transmission loss due to clutter (when clutter obstacle placed
% at the first path point next to the terminal) as function of spacing
% Equidistant segmentation is used.

clear all
close all


linestyle = {'b-', 'r--', 'k-.', 'g--', 'm-.', 'c--'};
%% Link parameters

% Flat terrain flag
flat = false;
% Frequency (GHz)
f = 10;

% Time percentage
p = 50;

% Tx-Rx distance (km)
dist = 10;


% Tx and Rx antenna heights
htg = [10, 15, 19.9, 20, 20.1];

hrg = htg;

% hci is the height of the first clutter point
hc1 = 20;   % Transmitter side (m)
hc2 = 20;   % Receiver side (m)

dstep = linspace(30,100,100);



% Additional input parameters for P.452-17
phi_path = 48;
Gt = 0;
Gr = 0;
pol = 1;
dct = 500;
dcr = 500;
DN = 45;
N0 = 350;
press = 1013;
temp = 20;


for hh = 1:length(htg)
    
    % Compute total transmission loss as a function of dc1=dc2
    for kk = 1:100
       
        nn = floor(dist/dstep(kk)*1000)+1;
        d = linspace(0,dist,nn);
        h = zeros(size(d));
        dc1 = d(2)-d(1);
        dc2 = d(2)-d(1);

        if (~flat)
            % include terrain - ideal triangular hill
            hmax = 1000;
            mid = floor(length(d)/2);
            mid1 = floor(mid/2);
            mid2 = floor(3*mid/2);
            for i = mid1 : mid
                h(i) = (i-mid1)*hmax/(mid-mid1);
            end
            for i = mid+1 : mid2
                h(i) = hmax - (i-mid)*hmax/(mid2-mid);
            end
        end

        % Initialize vector of clutter heights (no clutter)
        g = h;

        % Vector of zones (inland only)
        zone = 2*ones(size(d));

        
        dnew = d;
        gnew = g;

        % move the clutter within the first segment
        gnew(2) =  hc1+h(2);
        % symmetrically move the clutter within the last segement
        gnew(end-1) = hc2+h(end-1);



        % Compute basic transmission loss with clutter
        Lb_c(kk) = tl_p452_pdr(f, ...
            p, ...
            dnew, ...
            h, ...
            zone, ...
            gnew, ... % clutter + terrain profile along the path
            htg(hh), ...
            hrg(hh), ...
            phi_path,...
            Gt, ...
            Gr, ...
            pol, ...
            dct, ...
            dcr, ...
            DN, ...
            N0, ...
            press, ...
            temp);

        % Compute basic transmission loss without clutter
        Lb(kk) = tl_p452_pdr(f, ...
            p, ...
            dnew, ...
            h, ...
            zone, ...
            h, ... % 0 + terrain profile along the path
            htg(hh), ...
            hrg(hh), ...
            phi_path,...
            Gt, ...
            Gr, ...
            pol, ...
            dct, ...
            dcr, ...
            DN, ...
            N0, ...
            press, ...
            temp);


    end

        plot(dstep, Lb_c-Lb, linestyle{hh}, 'LineWidth', 2);
    hold on

    legendstr{hh}= [ 'h_{Tx} = h_{Rx} = ' num2str(htg(hh)) ];
    
end


if(flat)
    titlestr = ['Flat terrain,  h_{c1} = h_{c2} = ' num2str(hc1) ' m, d = ' num2str(dist) ' km'  ];
else
    titlestr = ['H = ' num2str(hmax) ' m, h_{c1} = h_{c2} = ' num2str(hc1) ' m, d = ' num2str(dist) ' km' ];
end



xlabel('d_{c1} = d_{c2} (m)')
ylabel('Excess L_b (dB)');

title(titlestr);

grid on

legend(legendstr, 'Location', 'Best')
set(gca, 'FontSize', 12)
figure

%% Transmission loss as a function of distance from Tx along the profile
% no clutter applied

for hh= 1:length(htg)
    for i = 10:length(d)

        dcrop = d(1:i);
        hcrop = h(1:i);
        zcrop = zone(1:i);

        dd(i-9) = dcrop(end);
        % Compute basic transmission loss without clutter
        Lbcrop(hh,i-9) = tl_p452_pdr(f, ...
            p, ...
            dcrop, ...
            hcrop, ...
            zcrop, ...
            hcrop, ... % 0 + terrain profile along the path
            htg(hh), ...
            hrg(hh), ...
            phi_path,...
            Gt, ...
            Gr, ...
            pol, ...
            dct, ...
            dcr, ...
            DN, ...
            N0, ...
            press, ...
            temp);



    end
    plot(dd, Lbcrop(hh,:), linestyle{hh},'LineWidth', 2);
    hold on
end

if(flat)
    titlestr = ['Flat terrain,  d = ' num2str(dist) ' km,' ' \Deltad = ' num2str(dstep(end)) ' m.' ];
else
    titlestr = ['H = ' num2str(hmax) ' m, d = ' num2str(dist) ' km,' ' \Deltad = ' num2str(dstep(end)) ' m.' ];
end


xlabel('d (km)')
ylabel('L_b (dB)')
legend(legendstr, 'Location', 'Best')
title(titlestr);
set(gca, 'FontSize', 12)
grid on

%% Transmission loss as a function of distance from Tx along the profile
% clutter applied (at a distance corresponding to the profile segment)

figure
for hh= 1:length(htg)
    for i = 10:length(d)

        dcrop = d(1:i);
        hcrop = h(1:i);
        zcrop = zone(1:i);
        gcrop = hcrop;
        gcrop(2) = hc1 + hcrop(2);
        gcrop(end-1) = hc2 + hcrop(end-1);


        dd(i-9) = dcrop(end);
        % Compute basic transmission loss without clutter
        Lbcrop(hh,i-9) = tl_p452_pdr(f, ...
            p, ...
            dcrop, ...
            hcrop, ...
            zcrop, ...
            gcrop, ... % 0 + terrain profile along the path
            htg(hh), ...
            hrg(hh), ...
            phi_path,...
            Gt, ...
            Gr, ...
            pol, ...
            dct, ...
            dcr, ...
            DN, ...
            N0, ...
            press, ...
            temp);



    end
    plot(dd, Lbcrop(hh,:), linestyle{hh},'LineWidth', 2);
    hold on
end

if(flat)
    titlestr = ['Flat terrain,  d = ' num2str(dist) ' km,' ' h_{c1} = ' num2str(hc1) ' m, h_{c2} = ' num2str(hc2) ' m, \Deltad = ' num2str(dstep(end)) ' m.' ];
else
    titlestr = ['H = ' num2str(hmax) ' m, d = ' num2str(dist) ' km, h_{c1} = ' num2str(hc1) ' m, h_{c2} = ' num2str(hc2) ' m, \Deltad = ' num2str(dstep(end)) ' m.' ];
end


xlabel('d (km)')
ylabel('L_{bc} (dB)')
legend(legendstr, 'Location', 'Best')
title(titlestr);
set(gca, 'FontSize', 12)
grid on


%% Plot the terrain profile
figure
plot(d,h,'b', 'LineWidth', 2);
hold on
plot(d,g,'r--','LineWidth', 2)
set(gca, 'FontSize', 12)
xlabel('d (km)')
ylabel('h (m)')
grid on
set(gca,'FontSize', 12)
