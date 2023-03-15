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

k1 = find( d<= d_tcv);
i1 = k1(end);

k2 = find(d >= d_tcv);
i2 = k2(1);

% apply linear interpolation
d1 = d(i1);
d2 = d(i2);
h1 = h(i1);
h2 = h(i2);
ds = d_tcv;

if (d1 == d2)
    hs = h1;
else
    hs = h1 + (h2-h1)*(ds-d1)/(d2-d1);
end

return
end
