function [Phipnte, Phipntn, Bt2r, dgc] = great_circle_path(Phire, Phite, Phirn, Phitn, Re, dpnt)
%great_circle_path Great-circle path calculations according to Attachment H
%   This function computes the great-circle intermediate points on the
%   radio path as defined in ITU-R P.2001-2 Attachment H
%
%     Input parameters:
%     Phire   -   Receiver longitude, positive to east (deg)
%     Phite   -   Transmitter longitude, positive to east (deg)
%     Phirn   -   Receiver latitude, positive to north (deg)
%     Phitn   -   Transmitter latitude, positive to north (deg)
%     Re      -   Average Earth radius (km)
%     dpnt    -   Distance from the transmitter to the intermediate point (km)
%
%     Output parameters:
%     Phipnte -   Longitude of the intermediate point (deg)
%     Phipntn -   Latitude of the intermediate point (deg)
%     Bt2r    -   Bearing of the great-circle path from Tx towards the Rx (deg)
%     dgc     -   Great-circle path length (km)
%
%     Example:
%     [Bt2r, Phipnte, Phipntn, dgc] = great_circle_path(Phire, Phite, Phirn, Phitn, Re, dpnt)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    12JUL16     Ivica Stevanovic, OFCOM         Initial version

%% H.2 Path length and bearing

% Difference (deg) in longitude between the terminals (H.2.1)

Dlon = Phire - Phite;

% Calculate quantity r (H.2.2)

r = sind(Phitn) * sind(Phirn) + cosd(Phitn) * cosd(Phirn) * cosd(Dlon);

% Calculate the path length as the angle subtended at the center of
% average-radius Earth (H.2.3)

Phid = acos(r);  % radians

% Calculate the great-circle path length (H.2.4)

dgc = Phid * Re;  % km

% Calculate the quantity x1 (H.2.5a)

x1 = sind(Phirn)-r*sind(Phitn);

% Calculate the quantity y1 (H.2.5b)

y1 = cosd(Phitn)*cosd(Phirn)*sind(Dlon);

% Calculate the bearing of the great-circle path for Tx to Rx (H.2.6)

if (abs(x1) < 1e-9 && abs(y1) < 1e-9 )
    Bt2r = Phire;
else
    Bt2r = atan2d(y1,x1);
end

%% H.3 Calculation of intermediate path point

% Calculate the distance to the point as the angle subtended at the center
% of average-radius Earth (H.3.1)

Phipnt = dpnt/Re;  %radians

% Calculate quantity s (H.3.2)

s = sind(Phitn)*cos(Phipnt) + cosd(Phitn)*sin(Phipnt)*cosd(Bt2r);

% The latitude of the intermediate point is now given by (H.3.3)

Phipntn = asind(s); % degs

% Calculate the quantity x2 (H.3.4a)

x2 = cos(Phipnt)-s*sind(Phitn);

% Calculate the quantity y2 (H.3.4b)

y2 = cosd(Phitn)*sin(Phipnt)*sind(Bt2r);

% Calculate the longitude of the intermediate point Phipnte (H.3.5)

if (x2 < 1e-9 && y2 < 1e-9)
    Phipnte = Bt2r;
else
    Phipnte = Phite + atan2d(y2,x2);
end

return
end
