function Ldft = dl_se_ft(d, hte, hre, adft, f, omega)
%dl_se_ft First-term part of spherical-Earth diffraction according to ITU-R P.452-16
%   This function computes the first-term part of Spherical-Earth diffraction
%   loss exceeded for p% time for antenna heights
%   as defined in Sec. 4.2.2.1 of the ITU-R P.452-16
%
%     Input parameters:
%     d       -   Great-circle path distance (km)
%     hte     -   Effective height of interfering antenna (m)
%     hre     -   Effective height of interfered-with antenna (m)
%     adft    -   effective Earth radius (km)
%     f       -   frequency (GHz)
%     omega   -   fraction of the path over sea
%
%     Output parameters:
%     Ldft   -   The first-term spherical-Earth diffraction loss not exceeded for p% time
%                Ldft(1) is for the horizontal polarization
%                Ldft(2) is for the vertical polarization
%
%     Example:
%     Ldft = dl_se_ft(d, hte, hre, adft, f, omega)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23DEC15     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    16DEC16     Ivica Stevanovic, OFCOM         corrected bug in dl_se_ft_inner function call where epsr and sigma order was interchanged


%% Body of function

% First-term part of the spherical-Earth diffraction loss over land

epsr = 22;
sigma = 0.003;

Ldft_land = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f);

% First-term part of the spherical-Earth diffraction loss over sea

epsr = 80;
sigma = 5;

Ldft_sea = dl_se_ft_inner(epsr, sigma, d, hte, hre, adft, f);


% First-term spherical diffraction loss 

Ldft = omega * Ldft_sea + (1-omega)*Ldft_land;      % Eq (29)

return
end
