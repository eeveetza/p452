function struct = read_C3_1_profile(filename)
% this function reads the .csv file in ASCII format
% as provided by Stephen Solomon
%   Tx Name, Tx Country
%   Tx Latitude, Tx Longitude, Tx Ant Gain
%   Tx Ground Elev, TX_AHAG
%   Rx Name, Rx Country
%   Rx Latitude, Rx Longitude, Rx Ant Gain
%   Rx Ground Elev, RX_AHAG
%   Freq MHz
%   Distance(i), Height(i), Latitude(i), Longitude(i), GridSpacing(i)  for i=1:n
%
% and returns a structure containing path profile (d, h, lat, lon, 
%
% Inputs:
% Variable    Unit  Type           Description
% filename    -     string         name of the terrain profile in .csv format
% 
% Outputs: 
% Variable    Unit  Type           Description
% struct      -     structure      structure containing the following fields
% .tx.the     deg   double         Tx latitude
% .tx.phi     deg   double         Tx longitude
% .tx.g       dBi   double         Tx antenna gain
% .tx.hasl    m     double         Tx ground elevation
% .tx.ahag    m     double         Tx antenna height above ground 
% .tx.name    m     string         Tx name
% .tx.country m     string         Tx country
% .rx.the     deg   double         Tx latitude
% .rx.phi     deg   double         Rx longitude
% .rx.g       dBi   double         Rx antenna gain
% .rx.hasl    m     double         Rx ground elevation
% .rx.ahag    m     double         Rx antenna height above ground 
% .rx.name    m     string         Rx name
% .rx.country m     string         Rx country
% .f          GHz   double         Frequency
% .d          km    double[]       Vector of terrain profile distances
% .the        deg   double[]       Vector of terrain profile latitudes
% .phi        deg   double[]       Vector of terrain profile longitudes
% .h          m     double[]       Vector of terrain profile elevations   

%   Example:
%   struct = read_C3_1_profile('C_3_1_profiles/3_1_0001.csv'); 

%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    02MAR20     Ivica Stevanovic, OFCOM         Initial version


fid = fopen(filename);

if (fid == -1)
    error([filename ' does not exist.']);
    return
end

try
    
    dummy = split(fgetl(fid),',');
    struct.tx.name = dummy{1};
    struct.tx.country = dummy{2};    
    
    
    dummy = split(fgetl(fid),',');
    struct.tx.the = str2double(dummy{1});
    struct.tx.phi = str2double(dummy{2});
    struct.tx.g = str2double(dummy{3});
    
   
    dummy = split(fgetl(fid),',');
    struct.tx.hasl = str2double(dummy{1});
    struct.tx.ahag = str2double(dummy{2});
    
       
    % receiver
    
    dummy = split(fgetl(fid),',');
    struct.rx.name = dummy{1};
    struct.rx.country = dummy{2};    
    
    dummy = split(fgetl(fid),',');
    struct.rx.the = str2double(dummy{1});
    struct.rx.phi = str2double(dummy{2});
    struct.rx.g = str2double(dummy{3});
    
    dummy = split(fgetl(fid),',');
    struct.rx.hasl = str2double(dummy{1});
    struct.rx.ahag = str2double(dummy{2});
    
    dummy = split(fgetl(fid),',');
    struct.f = str2double(dummy{1})/1000; %frequency in GHz
    
    count = 1;
    while (1)
       linestr = fgetl(fid);
       if (linestr == -1)
           break;
       end
       
       dummy =  split(linestr,',');
       
       struct.d(count) = str2double(dummy{1});
       struct.h(count) = str2double(dummy{2});
       struct.the(count) = str2double(dummy{3});
       struct.phi(count) = str2double(dummy{4});
       count = count + 1;
    end
    
    
catch
    error([filename ' does not have the expected format.']);
end

fclose(fid);

return