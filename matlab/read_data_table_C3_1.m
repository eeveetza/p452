function stano = read_data_table_C3_1(filename)

fprintf(1,'\nImporting data from: %s\n', filename);

 column = 'A:BV'; 
 sheet = 'Sheet1';
        
 C = readcell(filename, 'Sheet', 'Sheet1', "NumHeaderLines", 1);
 
 C(cellfun(@(x) all(ismissing(x)), C)) = {NaN};

 [nr, nc] = size(C);
 
 for i = 1 : nr
     stano{i}.number = C{i,1};
     stano{i}.fGHz   = C{i,2};
     stano{i}.dkm    = C{i,3};

     stano{i}.tx.name = C{i, 5};
     stano{i}.tx.lat  = C{i,7};
     stano{i}.tx.lon  = C{i,8};
     if stano{i}.tx.lon > 180
         stano{i}.tx.lon = stano{i}.tx.lon-360;
     end
     stano{i}.tx.hasl  = C{i,9};
     stano{i}.tx.ahag  = C{i,10};
     stano{i}.tx.g  = C{i,11};
     
     stano{i}.rx.name = C{i,12};
     stano{i}.rx.lat  = C{i,14};
     stano{i}.rx.lon  = C{i,15};
     if stano{i}.rx.lon > 180
         stano{i}.rx.lon = stano{i}.rx.lon-360;
     end
     stano{i}.rx.hasl  = C{i,16};
     stano{i}.rx.ahag  = C{i,17};
     stano{i}.rx.g  = C{i,18};
     
     stano{i}.tx.z = C{i,21};
     stano{i}.rx.z = C{i,22};
     
     stano{i}.N0 = C{i,25};
     stano{i}.DN = C{i,29};
     stano{i}.hm = C{i,38};
     stano{i}.dlt = C{i,40};
     stano{i}.dlr = C{i,44};
     stano{i}.t   = [0.001 0.002 0.003 0.005 0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 1 2 3 5 10 20 30 50 90 99];
     stano{i}.btl   = [C{i,53:74}];
     
 end
 
 return
%      
%  filename_out = 'table_C3_1.mat';
%  save(filename_out, 'stano');