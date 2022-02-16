function [ Ldp, Ld50 ] = dl_p( d, h, hts, hrs, hstd, hsrd, f, omega, p, b0, DN )
%dl_p Diffraction loss model not exceeded for p% of time according to P.452-16
%   function [Ldp, Ld50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
%
%   This function computes the diffraction loss not exceeded for p% of time
%   as defined in ITU-R P.452-16 (Section 4.5.4)
%
%     Input parameters:
%     d       -   vector of distances di of the i-th profile point (km)
%     h       -   vector of heights hi of the i-th profile point (meters
%                 above mean sea level. Both vectors contain n+1 profile points
%     hts     -   transmitter antenna height in meters above sea level (i=0)
%     hrs     -   receiver antenna height in meters above sea level (i=n)
%     hstd    -   Effective height of interfering antenna (m amsl) c.f. 5.1.6.3
%     hsrd    -   Effective height of interfered-with antenna (m amsl) c.f. 5.1.6.3
%     f       -   frequency expressed in GHz
%     omega   -   the fraction of the path over sea
%     p       -   percentage of time
%     b0      -   the time percentage that the refractivity gradient (DELTA-N) exceeds 100 N-units/km in the first 100m of the lower atmosphere
%     DN      -   the average radio-refractive index lapse-rate through the
%                 lowest 1 km of the atmosphere. Note that DN is positive
%                 quantity in this procedure
%
%     Output parameters:
%     Ldp    -   diffraction loss for the general path not exceeded for p % of the time 
%                according to Section 4.2.4 of ITU-R P.452-16. 
%                Ldp(1) is for the horizontal polarization 
%                Ldp(2) is for the vertical polarization
%     Ld50   -   diffraction loss for p = 50%
%
%     Example:
%     [Ldp, Ld50] = dl_p( d, h, hts, hrs, hstd, hsrd, ap, f, omega, p, b0, DN )
%       
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    01JAN16     Ivica Stevanovic, OFCOM         Initial version


%% Body of function

% Use the method in 4.2.3 to calculate diffractino loss Ld for effective 
% Earth radius ap = ae as given by equation (6a). Set median diffractino
% loss to Ldp50

[ae, ab] = earth_rad_eff(DN);

ap = ae;

Ld50 = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega );

if p == 50
    Ldp = Ld50;
    return
end

if p < 50
    
    % Use the method in 4.2.3 to calculate diffraction loss Ld for effective
    % Earth radius ap = abeta, as given in equation (6b). Set diffraction loss
    % not exceeded for beta0% time Ldb = Ld
    
    ap = ab;
    
    Ldb = dl_delta_bull( d, h, hts, hrs, hstd, hsrd, ap, f, omega );

    % Compute the interpolation factor Fi
    
    if p > b0
        
        Fi = inv_cum_norm(p/100) / inv_cum_norm(b0/100);   % eq (41a)
        
    else
        
        Fi = 1;
        
    end
    
    % The diffraction loss Ldp not exceeded for p% of time is now given by
    
    Ldp = Ld50 + Fi*(Ldb - Ld50);   % eq (42)
    
end

return
end