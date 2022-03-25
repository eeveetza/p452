function c = readcsv(filename)
%%
%     This function reads a csv file from filename and returns a cell array
%     it assumes that every
%
%     Input arguments:
%     filename  -   name of the file to read with comma separated values
%
%     Output arguments:
%     c   -   cell array containing the comma separated values from
%
%     Example:
%     c = readcsv(filename)
%
%     Rev   Date        Author                          Description
%     -------------------------------------------------------------------------------
%     v0    16FEB22     Ivica Stevanovic, OFCOM         First implementation in matlab

fid = fopen(filename, 'r');

if (fid==-1)
    error(['File ' filename ' cannot be found.']);
end
 
nrows = 0;
ncols = 1;

while (true)
    line = fgetl(fid);
    
    if (~ischar(line))
        break
    end
    
    nrows = nrows + 1;
    
    if nrows == 2
        %ncols = length(split(line,','));
         ncols = length( regexp(line,',','split') );
        end
end

c = cell(nrows-1,ncols);


fseek(fid,0,'bof');

for i = 1:nrows
    
   dummy = fgetl(fid);
   
   if (i == 1)  % do not read the header
       continue
   end
   %c(i-1,:) = split(dummy,',');
   c(i-1,:) = regexp(dummy,',','split');
end

fclose(fid);


return

end