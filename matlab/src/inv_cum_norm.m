function I = inv_cum_norm( x )
%inv_cum_norm approximation to the inverse cummulative normal distribution
%   I = inv_cum_norm( x )
%   This function implements an approximation to the inverse cummulative
%   normal distribution function for x <= 0.5 as defined in Attachment 3 to
%   Annex 1 of the ITU-R P.452-16
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02JAN16     Ivica Stevanovic, OFCOM         Initial version

if x < 0.000001
    x = 0.000001;
end

if x > 0.5
    warning('This function is defined for arguments not larger than 0.5');
end

tx = sqrt(-2*log(x));    % eq (172a)

C0 = 2.515516698;        % eq (172c)
C1 = 0.802853;           % eq (172d)
C2 = 0.010328;           % eq (172e)
D1 = 1.432788;           % eq (172f)
D2 = 0.189269;           % eq (172g)
D3 = 0.001308;           % eq (172h)

ksi = ( (C2*tx+C1)*tx + C0 )/ ( ((D3*tx + D2)*tx + D1)*tx + 1 );  % eq (172b)

I = ksi - tx;            % eq (172)

return 
end

