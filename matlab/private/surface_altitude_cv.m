function   hs = surface_altitude_cv(h, d, d_tcv)
%surface_altitude_cv altitude on the surface of the Earth below common volume
%     hs = surface_altitude_cv(h, d, d_tcv)
%     This function computes the altitude of the point at the surface of
%     the Earth below common volume
%
%     Input arguments:
%     d       -   vector of distances in the path profile (km)
%     h       -   vector of heights (masl)
%     d_ctv   -   horizontal path length from transmitter to common volume computed using (3.9.1a)
%
%     Output arguments:
%     hs      -   altitude on the surface of the Earth below common volume (masl)
%
%     Example:
%     hs = surface_altitude_cv(h, d, d_tcv)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    17FEB23     Ivica Stevanovic, OFCOM         First implementation in matlab

if (d_tcv > d(end) || d_tcv < 0)
    error('Horizontal distance of the common volume to the transmitter must be within  [0, d(end)]');
end

[~,ii] = min(abs(d-d_tcv));

if(d(ii) == d_tcv)
    i1 = ii;
    hs = h(i1);
    return
end

if d(ii) < d_tcv
    i1 = ii;
    i2 = ii + 1;
else
    i2 = ii;
    i1 = ii - 1;
end

% apply linear interpolation

    hs = h(i1) + (h(i2)-h(i1))*(d_tcv-d(i1))/(d(i2)-d(i1));

return
end
