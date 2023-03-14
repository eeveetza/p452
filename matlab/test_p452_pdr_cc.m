%% MATLAB script used to analyze the behavior of PDR P.452-clutter
% For selected corner cases
% Scenario 1: Flat terrain
% flat = true
% dist = 10, 30, 60 m (free-space, spherical-earth diffraction, tropospheric and other effects)
% Scenario 2: Ideal Triangular "Mountain"
% flat = false
% dist = 10, 30, 60
% The script computes and plots
% Transmission loss as a function of distance (w/o clutter)
% Excess transmission loss due to clutter (when clutter obstacle placed
% within the first segment)
% Equidistant segmentation is used except for the first two/ last two
% segments with terminal clutter.

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

% Step in terrain profile (m) % needs to be integer in this implementation
dstep = 40;

% Tx and Rx antenna heights
htg = [10, 15, 19.9, 20, 20.1];
hrg = htg;

% hci is the height of the first clutter point
hc1 = 20;   % Transmitter side (m)
hc2 = 20;   % Receiver side (m)

% dci is the distance of the first clutter point from the terminal
% at 1 m increment
dcmin = 1;
dc1 = linspace(dcmin,dstep,dstep)/1000;
dc2 = dc1;

% Number of points in the profile
nstep = floor((dist-2*dstep/1000)/dstep*1000) + 1;

% Vector of distances in the path profile
d1 = union(0, dc1);
d2 = d1 + dist - dstep/1000;
d3 = linspace(dstep/1000, dist-dstep/1000, nstep);
d = union(d1,d3);
d = union(d, d2);
% Vector of heights in the path profile (assumed flat terrain 0 m)
h = zeros(size(d));

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
    for i = 2:length(dc1)+1

        dnew = d;
        gnew = g;

        % move the clutter within the first segment
        gnew(i) =  hc1+h(i);
        % symmetrically move the clutter within the last segement
        gnew(end-i+1) = hc2+h(end-i+1);

        % Compute basic transmission loss with clutter
        Lb_c(hh,i-1) = tl_p452_pdr(f, ...
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
        Lb(hh,i-1) = tl_p452_pdr(f, ...
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

    legendstr{hh}= [ 'h_{Tx} = h_{Rx} = ' num2str(htg(hh)) ];
    
end


if(flat)
    titlestr = ['Flat terrain,  h_{c1} = h_{c2} = ' num2str(hc1) ' m, d = ' num2str(dist) ' km,' ' \Deltad = ' num2str(dstep) ' m.' ];
else
    titlestr = ['H = ' num2str(hmax) ' m, h_{c1} = h_{c2} = ' num2str(hc1) ' m, d = ' num2str(dist) ' km,' ' \Deltad = ' num2str(dstep) ' m.' ];
end


figure
for (hh = 1:length(htg))
    plot(dc1*1000, Lb_c(hh,:)-Lb(hh,:), linestyle{hh}, 'LineWidth', 2);
    hold on
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
    titlestr = ['Flat terrain,  d = ' num2str(dist) ' km,' ' \Deltad = ' num2str(dstep) ' m.' ];
else
    titlestr = ['H = ' num2str(hmax) ' m, d = ' num2str(dist) ' km,' ' \Deltad = ' num2str(dstep) ' m.' ];
end

xlabel('d (km)')
ylabel('L_b (dB)')
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
