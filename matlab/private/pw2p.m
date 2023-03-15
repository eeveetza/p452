function   p = pw2p(pw, phi, db, d)
%pw2p Annual equivalent time percentage from the worst-month time percentage
%     p = pw2p(pw, phi, db, d)
%     This function computes the annual equivalent time percentage p
%     starting from the worst-month time percentage, pw
%     as defined in ITU-R P.452-16.
%
%     Input arguments:
%     pw      -   the worst-month time percentage
%     phi     -   path centre latitude (deg)
%     db      -   aggregate length of the path sections over water (km)
%     d       -   Great-circle path distance (km) calculated using (148)
%
%     Output arguments:
%     p       -   annual equivalent time percentage
%
%     Example:
%     p = pw2p(pw, phi, db, d)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    16MAR23     Ivica Stevanovic, OFCOM         Corrected bugs as indicated by mplussier
%                                                       Added capping from below set to pw/12    


% Fraction of the path over water

omega = db/d;

if abs(phi) <= 45
   GL = sqrt(1.1+(abs(cosd(2*phi))).^0.7);    %(1a)
else
   GL = sqrt(1.1-(abs(cosd(2*phi))).^0.7);
end

pexp = ( log10(pw) + log10(GL)-0.186*omega - 0.444 ) / ...    %(1)
        (0.816 + 0.078*omega);

p = 10.^pexp;

% If necessary the value of p must be limited such that p >= pw/12
p = max(p, pw/12);

return
end