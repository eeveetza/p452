function stano = read_data_table_C3_1(filename)

fprintf(1,'\nImporting data from: %s\n', filename);

 column = 'A:BV'; 
 sheet = 'Sheet1';
        
 [num, txt, raw] =  xlsread(filename, sheet, column);
 
 [nr, nc] = size(num);
 
 for i = 1 : nr
     stano{i}.number = num(i,1);
     stano{i}.fGHz   = num(i,2);
     stano{i}.dkm    = num(i,3);

     stano{i}.tx.name = string(txt(i+1,5));
     stano{i}.tx.lat  = num(i,7);
     stano{i}.tx.lon  = num(i,8);
     if stano{i}.tx.lon > 180
         stano{i}.tx.lon = stano{i}.tx.lon-360;
     end
     stano{i}.tx.hasl  = num(i,9);
     stano{i}.tx.ahag  = num(i,10);
     stano{i}.tx.g  = num(i,11);
     
     stano{i}.rx.name = string(txt(i+1,12));
     stano{i}.rx.lat  = num(i,14);
     stano{i}.rx.lon  = num(i,15);
     if stano{i}.rx.lon > 180
         stano{i}.rx.lon = stano{i}.rx.lon-360;
     end
     stano{i}.rx.hasl  = num(i,16);
     stano{i}.rx.ahag  = num(i,17);
     stano{i}.rx.g  = num(i,18);
     
     stano{i}.tx.z = string(txt(i+1,21));
     stano{i}.rx.z = string(txt(i+1,22));
     
     stano{i}.N0 = num(i,25);
     stano{i}.DN = num(i,29);
     stano{i}.hm = num(i,38);
     stano{i}.dlt = num(i,40);
     stano{i}.dlr = num(i,44);
     stano{i}.t   = [0.001 0.002 0.003 0.005 0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 90 99];
     stano{i}.btl   = num(i,53:74);
     
 end
 
 return
%      
%  filename_out = 'table_C3_1.mat';
%  save(filename_out, 'stano');