function   omega = path_fraction(d, zone, zone_r)
%path_fraction Path fraction belonging to a given zone_r
%     omega = path_fraction(d, zone, zone_r)
%     This function computes the path fraction belonging to a given zone_r
%     of the great-circle path (km) 
%
%     Input arguments:
%     d       -   vector of distances in the path profile
%     zone    -   vector of zones in the path profile
%     zone_r  -   reference zone for which the fraction is computed
%
%     Output arguments:
%     omega   -   path fraction belonging to the given zone_r
%
%     Example:
%     omega = path_fraction(d, zone, zone_r)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02FEB16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    10NOV22     Ivica Stevanovic, OFCOM         Corrected a bug in line 36 (suggested by Martin-Pierre Lussier @mplussier)

dm = 0;

[start,stop] = find_intervals((zone == zone_r));

n = length(start);

for i = 1:n
    delta = 0;
    if (d(stop(i))<d(end))
        delta = delta + ( d(stop(i)+1)-d(stop(i)) )/2.0;
    end
    
    if (d(start(i))>0)
        delta = delta + ( d(start(i))-d(start(i)-1) )/2.0;
    end
    
   dm = dm + d(stop(i))-d(start(i)) + delta;
   
end

omega = dm/(d(end)-d(1));

return
end