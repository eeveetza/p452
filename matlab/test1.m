d = linspace(0,1000,1001);
h = zeros(size(d));
z = ones(size(d));
GHz = 10;
Tpc = 50;
Phire = 70.5;
Phirn = -60;
Phite = 75;
Phitn = -60;
Re = 6371;


Hrg = 100;
Htg = 100;
tropo = false; %when set to true, returns only troposcatter basic tl

Grx = 0;
Gtx = 0;
pol = 1;
dct = 500;
dcr = 500;
press = 1013;
temp = 20;

count = 1;
imin = 6;
for ii = imin:length(d)
    dd = d(1:ii);
    hh = h(1:ii);
    zz = z(1:ii);

    L0(count) = tl_p452_pdr(GHz, Tpc, dd, hh, hh, zz, Htg, Hrg, Phite, Phitn, Phire, Phirn, Gtx, Grx, pol, dct, dcr, press, temp, 0, false);
    L1(count) = tl_p452_pdr(GHz, Tpc, dd, hh, hh, zz, Htg, Hrg, Phite, Phitn, Phire, Phirn, Gtx, Grx, pol, dct, dcr, press, temp, 1, false);
    Lfs(count) = 92.4 + 20*log10(GHz) + 10*log10(dd(end).^2 + (Htg-Hrg).^2/1e6);
    count = count + 1;
end

plot(d(imin:end), L0, 'b', 'LineWidth', 2)
hold on
plot(d(imin:end), L1, 'r', 'LineWidth', 2)
grid on
plot(d(imin:end), Lfs, 'g', 'LineWidth', 2)
set(gca,'FontSize', 14)
legend('P.452-18','PDR P.452 ', 'Free-Space')
xlabel('distance (km)')
ylabel('Lb (dB)')
titlestr = ['f = ' num2str(GHz) ' GHz, Htg = ' num2str(Htg) ' m, Hrg = ' num2str(Hrg) ' m' ]
title(titlestr)