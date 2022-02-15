function [lat, lon] = gcintermediate(lat1, lon1, lat2, lon2, f)
%gcintermediate computes an intermediate point on a great circle
%     [lat, lon] = gcintermediate(lat1, lon1, lat2, lon2, f)
%
%     This function computes an intermediate point on a great circle
%     on a path between the points with coordinates (lat1, lon1) and (lat2, lon2)
%     and for a given fraction f of distance between those two points
%     f=0 defines point 1. f=1 defines point 2.
%
%     Input arguments:
%     lat1    -   latitude of the point 1 (deg)
%     lon1    -   longitude of the point 1 (deg)
%     lat2    -   latitude of the point 2 (deg)
%     lon2    -   longitude of the point 2 (deg)
%     f       -   fraction of the distance between the points 1 and 2
%
%     Output arguments:
%     lat     -   latitude of the intermediate point (deg)
%     lon     -   longitude of the intermediate point (deg)
%
%     Example:
%     [lat, lon] = gcintermediate(lat1, lon1, lat2, lon2, f)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    23AUG16     Ivica Stevanovic, OFCOM         Initial version


d = 2*asind(sqrt((sind((lat1-lat2)/2))^2 + ...
    cosd(lat1)*cosd(lat2)*(sind((lon1-lon2)/2))^2));


A = sind((1-f)*d)/sind(d);
B = sind(f*d)/sind(d);
x = A*cosd(lat1)*cosd(lon1) +  B*cosd(lat2)*cosd(lon2);
y = A*cosd(lat1)*sind(lon1) +  B*cosd(lat2)*sind(lon2);
z = A*sind(lat1)            +  B*sind(lat2);
lat = atan2d(z,sqrt(x^2+y^2));
lon = atan2d(y,x);

             
return
end