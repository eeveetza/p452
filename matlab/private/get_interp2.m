function y = get_interp2(mapstr, phie, phin)
%get_inter2 Interpolates the value from Map at a given phie,phin
%
%     Input parameters:
%     map     -   string pointing to the radiometeorological map
%                 'DN50', or 'N050'
%                 (rows-latitude: 90 to -90, columns-longitude: 0 to 360)
%     phie    -   Longitude, positive to east (deg)
%     phin    -   Latitude, positive to north (deg)
%     spacing -   Resolution in latitude/longitude (deg)
%
%     Output parameters:
%     y      -    Interpolated value from the radiometeorological mapt at the point (phie,phin)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    17APR23     Ivica Stevanovic, OFCOM         Initial version

switch mapstr
    case 'DN50'
        map = DigitalMaps_DN50();
        nr = 121;
        nc = 241;
        spacing = 1.5;
    case 'N050'
        map = DigitalMaps_N050();
        nr = 121;
        nc = 241;
        spacing = 1.5;
end

longitudeOffset = phie;

if (phie < 0.0)
    longitudeOffset = phie + 360;
end

latitudeOffset = 90.0 - phin;

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
