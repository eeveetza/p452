function [d, h, zone] = test_profile(profile)

d = [];
h = [];
zone = [];

% Zone type: Coastal land; Code: A1; Code in this routine: 1
% Zone type: Inland      ; Code: A2; Code in this routine: 2
% Zone type: Sea         ; Code:  B; Code in this routine: 3

if (profile == 1)
    
    % The path profile as defined in the Excel worksheet Rec ITU-R P_452-14 (WP 3M Nov 2010)
    
    d = linspace(0,109,110);
    
    h = [40.000 24.000 35.000 38.000 41.000 39.000 41.000 41.000 44.000 35.000 ...
    39.000 24.000 26.000 32.000 36.000 43.000 49.000 48.000 48.000 45.000 41.000 ...
    35.000 36.000 36.000 42.000 43.000 52.000 55.000 73.000 67.000 68.000 53.000 ...
    50.000 14.000 1.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 ...
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 ...
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 ...
    0.000 0.000 0.000 0.000 0.000 0.000 0.000 45.000 54.000 63.000 79.000 85.000 39.000 ...
    21.000 15.000 37.000 42.000 63.000 51.000 35.000 68.000 68.000 53.000 63.000 ...
    75.000 92.000 174.000 184.000 145.000 158.000 122.000 129.000 97.000 100.000 ...
    91.000 144.000 119.000 154.000 183.000 ]; 
    
    for i = 1:35
        zone(i)   = 1;
    end
    for i = 36:78
        zone(i)   = 3;
    end
    for i = 79:97
        zone(i)   = 1;
    end
    for i = 98:103
        zone(i)   = 2;
    end
    zone(104)     = 1; 
    zone(105)     = 2; 
    zone(106)     = 1; 
    for i = 107:110
        zone(i)   = 2;
    end

elseif(profile==2)
    
    aa = load('profile2.dat');
    
    d = aa(:,1);
    h = aa(:,2);
    for i=1:length(d)
        zone(i) = 2;
    end
        
    
    elseif(profile==3)
    
    aa = load('profile3.dat');
    
    d = aa(:,1);
    h = aa(:,2);
    for i=1:length(d)
        zone(i) = 2;
    end
    
     elseif(profile==4)
    
    aa = load('profile4.dat');
    
    d = aa(:,1);
    h = aa(:,2);
    for i=1:length(d)
        zone(i) = 2;
    end
    
    elseif(profile==5)
    
    aa = load('profile5.dat');
    
    d = aa(:,1);
    h = aa(:,2);
    for i=1:length(d)
        zone(i) = 2;
    end
    
    
     elseif(profile==6)
    
    aa = load('profile6.dat');
    
    d = aa(:,1);
    h = aa(:,2);
    for i=1:length(d)
        zone(i) = 2;
    end
end

return
end