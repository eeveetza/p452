function climzone = cz_p1812_to_p452(climzone1)
%% cz_p1812_to_p452 transposes climatic zone codes of P1812 to P452
%
%           coastal     inland      sea
% p1812     3           4           1
% p452      1           2           3

climzone = climzone1;
for i = 1:length(climzone)

    switch climzone1(i)
        case 3
            climzone(i) = 1;
        case 4
            climzone(i) = 2;
        case 1
            climzone(i) = 3;
        otherwise
            error('Unknown climatic zone in function argument %d', climzone1);
    end
end
return
end