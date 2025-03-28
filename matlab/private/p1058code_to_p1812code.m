function p1812code = p1058code_to_p1812code(cc, ch)
%% p1058code_to_p1812code the coverage code transformation from P1058 to P1812
% This maping transforms the coverage code in P.1058 to the one in P.1812.

% 1 - Water/sea, 2 - Open/rural, 3 - Suburban, 4 - Urban/trees/forest, 5 - Dense urban
% 0 m          , 0 m           , 10 m        , 15 m                  , 20 m
% 60, 70       , 10, 11, 37    , 32, 33, 34, , 19, 20, 35, 37        , 36

p1812code = zeros(size(cc));
for i = 1:length(cc)
    switch cc(i)
        case 60
            p1812code(i) = 1;
        case 70
            p1812code(i) = 1;
        case 10
            p1812code(i) = 2;
        case 11
            p1812code(i) = 2;
        case 12
            p1812code(i) = 2;
        case 13
            p1812code(i) = 2;
        case 37
            if ch(i) == 0
                p1812code(i) = 2;
            else
                p1812code(i) = 4;
            end
        case 32
            p1812code(i) = 3;
        case 33
            p1812code(i) = 3;
        case 34
            p1812code(i) = 3;
        case 19
            p1812code(i) = 4;
        case 20
            p1812code(i) = 4;
        case 35
            p1812code(i) = 4;
        case 36
            p1812code(i) = 5;
        otherwise
            error(['Unknown clutter code ' num2str(cc(i))]);
    end


end

return
end
