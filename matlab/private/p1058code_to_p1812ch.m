function R = p1058code_to_p1812ch(cc, ch)
%% p1058code_to_p1812ch the coverage code transformation from P1058 to P1812 representetive clutter heights
% This maping transforms the coverage code in P.1058 to the representative clutter heights in P.1812.

% 1 - Water/sea, 2 - Open/rural, 3 - Suburban, 4 - Urban/trees/forest, 5 - Dense urban
% 0 m          , 0 m           , 10 m        , 15 m                  , 20 m
% 60, 70       , 10, 11, 12, 13    , 32, 33, 34, , 19, 20, 35, 37        , 36

R = zeros(size(cc));
for i = 1:length(cc)
    switch cc(i)
        case 60
            R(i) = 0;
        case 70
            R(i) = 0;
        case 10
            R(i) = 0;
        case 11
            R(i) = 0;
        case 12
            R(i) = 0;
        case 13
            R(i) = 0;
    
        case 32
            R(i) = 10;
        case 33
            R(i) = 10;
        case 34
            R(i) = 10;
        case 19
            R(i) = 15;
        case 20
            R(i) = 15;
        case 35
            R(i) = 15;
        case 37
            R(i) = 15;
        case 36
            R(i) = 20;
        otherwise
            error(['Unknown clutter code ' num2str(cc(i))]);
    end


end

return
end
