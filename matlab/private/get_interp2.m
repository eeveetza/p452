function y = get_interp2(mapstr, phie, phin)
%get_inter2 Interpolates the value from Map at a given phie,phin
%
%     Input parameters:
%     mapstr  -   string pointing to the radiometeorological map
%     phie    -   Longitude, positive to east (deg) [-180, 180]
%     phin    -   Latitude, positive to north (deg) [-90, 90]
%     spacing -   Resolution in latitude/longitude (deg)
%
%     Output parameters:
%     y      -    Interpolated value from the radiometeorological mapt at the point (phie,phin)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    09SEP24     Ivica Stevanovic, OFCOM         Initial version


if (phin < -90 || phin > 90)
    error ('Latitude must be within the range -90 to 90 degrees');
end

if (phie < -180 || phie > 180)
    error('Longitude must be within the range -180 to 180 degrees');
end

errorstr = sprintf(['DigitalMaps_%s() not found. \n' ...
    '\nBefore running get_interp2, make sure to: \n' ...
    '    1. Download and extract the required maps to ./private/maps:\n' ...
    '        - From ITU-R P.452-18: DN50.TXT and N050.TXT\n' ...
    '    2. Run the script initiate_digital_maps.m to generate the necessary functions.\n'], mapstr);


switch mapstr
    case 'DN50'
        try
            map = DigitalMaps_DN50();
        catch
            error(errorstr);
        end
        nr = 121;
        nc = 241;
        spacing = 1.5;
        % lat starts with 90
        latitudeOffset = 90 - phin;
        % lon starts with 0
        longitudeOffset = phie;
        if phie < 0
            longitudeOffset = phie + 360;
        end


    case 'N050'
        try
            map = DigitalMaps_N050();
        catch
            error(errorstr);
        end
        nr = 121;
        nc = 241;
        spacing = 1.5;
        % lat starts with 90
        latitudeOffset = 90 - phin;
        % lon starts with 0
        longitudeOffset = phie;
        if phie < 0
            longitudeOffset = phie + 360;
        end

    otherwise

        error('Error in function call. Uknown map: %s.\n',mapstr);
end

latitudeIndex  = floor(latitudeOffset / spacing)  + 1;
longitudeIndex = floor(longitudeOffset / spacing) + 1;

latitudeFraction  = (latitudeOffset / spacing)  - (latitudeIndex  - 1);
longitudeFraction = (longitudeOffset / spacing) - (longitudeIndex - 1);

val_ul = map(latitudeIndex, longitudeIndex);
val_ur = map(latitudeIndex, min(longitudeIndex + 1, nc));
val_ll = map(min(latitudeIndex + 1, nr), longitudeIndex);
val_lr = map(min(latitudeIndex + 1, nr), min(longitudeIndex + 1, nc));

y1 = longitudeFraction  * ( val_ur - val_ul ) + val_ul;
y2 = longitudeFraction  * ( val_lr - val_ll ) + val_ll;
y  = latitudeFraction * ( y2 - y1 ) + y1;

return
end
