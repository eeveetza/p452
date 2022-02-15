function [k1, k2] = find_intervals(series)
%find_intervals Find all intervals with consecutive 1's
%     [k1, k2] = find_intervals(series)
%     This function finds all 1's intervals, namely, the indices when the
%     intervals start and where they end
%
%     For example, for the input indices
%           0 0 1 1 1 1 0 0 0 1 1 0 0
%     this function will give back
%       k1 = 3, 10
%       k2 = 6, 11
%
%     Input arguments:
%     indices -   vector containing zeros and ones
%
%     Output arguments:
%     k1      -   vector of start-indices of the found intervals
%     k2      -   vector of end-indices of the found intervals
%
%     Example:
%     [k1, k2] = find_intervals(indices)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    22JAN16     Ivica Stevanovic, OFCOM         First implementation in matlab
%     v1    22JUN16     Roger LeClair, leclairtelecom   Modified to optimize for speed

k1 = [];
k2 = [];

if max(series) == 1
    k1 = find(diff([0, series]) == 1);
    k2 = find(diff([series, 0]) == -1);
    
end % if

return
end

% >>> Start code change.

%count = 1;

%if (series(1) ==1 && series(2) == 1)
%    k1(count) = 1;
%    count = count + 1;
%end

%for i = 2:length(series)-1
    
    
%    if (series(i-1) == 0 ) && (series(i) == 1)  
%        k1(count) = i;
%        count = count + 1;
        
%    end
%end

%count = 1;

%for i = 2:length(series)-1
   
%    if (series(i) == 1) && (series(i+1) == 0 )
%        k2(count) = i;
%        count = count + 1;
        
%    end
%end

%if (series(end) ==1 && series(end-1) == 1)
%    k2(count) = length(series);
%    count = count + 1;
%end

% >>> End code change.


