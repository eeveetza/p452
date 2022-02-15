function   b0 = beta0(phi, dtm, dlm)
%%
%     This function computes the time percentage for which refractive index
%     lapse-rates exceeding 100 N-units/km can be expected in the first 100
%     m of the lower atmosphere
%     as defined in ITU-R P.452-16.
%
%     Input arguments:
%     phi     -   path centre latitude (deg)
%     dtm     -   the longest continuous land (inland + coastal) section of the great-circle path (km)
%     dlm     -   the longest continuous inland section of the great-circle path (km)
%
%     Output arguments:
%     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
%
%     Example:
%     b0 = beta0(phi, dtm, dlm)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab


tau = 1- exp(-(4.12*1e-4*dlm.^2.41));       % (3a)

mu1 = ( ...                                 % (3)
       10^(-dtm/(16-6.6*tau)) + ...
       10^(-5*(0.496 + 0.354*tau)) ...
       ).^0.2;

if mu1 > 1
    mu1 = 1;
end
   
if abs(phi) <= 70
   mu4 = 10^( (-0.935 + 0.0176*abs(phi))*log10(mu1) );   % (4)
   b0 = 10^( -0.015*abs(phi) + 1.67 )*mu1*mu4;           % (2)   
else
   mu4 = 10^(0.3*log10(mu1));                            % (4)
   b0 = 4.17*mu1*mu4;                                    % (2)
end

return
end
