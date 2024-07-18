latcnt = 90:-1.5:-90;               %Table 2.4.1
loncnt = 0:1.5:360;                 %Table 2.4.1
[LON,LAT] = meshgrid(loncnt, latcnt);


Phin = (2*rand(180,1)-1)*90;
Phie = rand(360,1)*360;

dN1 = zeros(length(Phin),length(Phie));
dN2 = zeros(length(Phin),length(Phie));

N1 = zeros(length(Phin),length(Phie));
N2 = zeros(length(Phin),length(Phie));


        DN50 = DigitalMaps_DN50();
        N050 = DigitalMaps_N050();


tic
for nn = 1:length(Phin)
    for ee = 1:length(Phie)
        phim_e = Phie(ee);
        phim_n = Phin(nn);


        % Map phicve (-180, 180) to loncnt (0,360);
        phim_e1 = phim_e;
        if phim_e1 < 0
            phim_e1 = phim_e + 360;
        end

        dN1(nn,ee) = interp2(LON,LAT,DN50,phim_e1,phim_n);
        N1(nn,ee)  = interp2(LON,LAT,N050,phim_e1,phim_n);
        
    end
end
toc

tic
for nn = 1:length(Phin)
    for ee = 1:length(Phie)
        phim_e = Phie(ee);
        phim_n = Phin(nn);
        dN2(nn,ee) = get_interp2('DN50', phim_e,phim_n);
        N2(nn,ee)  = get_interp2('N050', phim_e,phim_n);
    end
end
toc

max(max(abs(N2-N1)))
max(max(abs(dN2-dN1)))