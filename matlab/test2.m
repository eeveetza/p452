% testing against ITM
d = union(linspace(0,100,101), linspace(100, 1000, 91));
h = zeros(size(d));
z = ones(size(d));
GHz = 10;
Tpc = 50;
Phire = 70.5;
Phirn = -60;
Phite = 75;
Phitn = -60;
Re = 6371;

Hrg = 1.5;
Htg = 1.5;
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

dpnt = dd(end)-dd(1);
[Phire_ii, Phirn_ii, bt2r, dgc] = great_circle_path(Phire, Phite, Phirn, Phitn, Re, dpnt);

    L0(count) = tl_p452(GHz, Tpc, dd, hh, zz, Htg, Hrg, Phite, Phitn, Phire_ii, Phirn_ii, Gtx, Grx, pol, dct, dcr, press, temp, 0);
    L1(count) = tl_p452(GHz, Tpc, dd, hh, zz, Htg, Hrg, Phite, Phitn, Phire_ii, Phirn_ii, Gtx, Grx, pol, dct, dcr, press, temp, 1);
    Lfs(count) = 92.4 + 20*log10(GHz) + 10*log10(dd(end).^2 + (Htg-Hrg).^2/1e6);
    count = count + 1;
end

plot(d(imin:end), L0, 'b', 'LineWidth', 2)
hold on
plot(d(imin:end), L1, 'r', 'LineWidth', 2)
grid on
plot(d(imin:end), Lfs, 'g', 'LineWidth', 2)
set(gca,'FontSize', 14)
legend('P.452-17','P.452 (new troposcatter)', 'Free-Space')
xlabel('distance (km)')
ylabel('Lb (dB)')
titlestr = ['f = ' num2str(GHz) ' GHz, Htg = ' num2str(Htg) ' m, Hrg = ' num2str(Hrg) ' m' ]
title(titlestr)