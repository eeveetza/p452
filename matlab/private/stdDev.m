function sigmaLoc = stdDev(f, h, R, wa)
%stdDev computes standard deviation of location variability
%   sigmaLoc = stdDev(f, h, R, wa)
%   This function computes the standard deviation according to ITU-R
%   P.1812-6, Annex 1, §4.8 and §4.10
%     Input parameters:
%     f       -   Frequency (GHz)
%     h       -   receiver/mobile antenna height above the ground (m)
%     R       -   height of representative clutter at the receiver/mobile location (m)
%     wa      -   prediction resolution, i.e., the width of the square area over which the variability applies

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    19MAR20     Ivica Stevanovic, OFCOM         Initial version

sigmaLoc = (0.52 + 0.024 * f) * wa.^0.28;

if (h < R)
    uh = 1;
else
    if (h >= R + 10)
        uh = 0;
    else
        uh = 1 - (h-R)/10;
    end
end

sigmaLoc = sigmaLoc * uh;

return
end



