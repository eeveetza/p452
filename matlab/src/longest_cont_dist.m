function   dm = longest_cont_dist(d, zone, zone_r)
%longest_cont_dist Longest continuous path belonging to the zone_r
%     dm = longest_cont_dist(d, zone, zone_r)
%     This function computes the longest continuous section of the
%     great-circle path (km) for a given zone_r
%
%     Input arguments:
%     d       -   vector of distances in the path profile
%     zone    -   vector of zones in the path profile
%     zone_r  -   reference zone for which the longest continuous section
%                 is computed
%
%     Output arguments:
%     dm      -   the longest continuous section of the great-circle path (km) for a given zone_r
%
%     Example:
%     dm = longest_cont_dist(d, zone, zone_r)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    12FEB16     Ivica Stevanovic, OFCOM         included zone_r==12
dm = 0;

if zone_r  == 12
    [start,stop] = find_intervals((zone == 1)+(zone==2));
else
    [start,stop] = find_intervals((zone == zone_r));
end

n = length(start);

for i = 1:n
    delta = 0;
    if (d(stop(i))<d(end))
        delta = delta + ( d(stop(i)+1)-d(stop(i)) )/2.0;
    end
    
    if (d(start(i))>0)
        delta = delta + ( d(stop(i))-d(stop(i)-1) )/2.0;
    end
    
   dm = max(d(stop(i))-d(start(i)) + delta, dm);
end


return
end