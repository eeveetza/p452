d = linspace(0,10000,1000);
h = zeros(size(d));
z = ones(size(d));
GHz = 10;
Tpc = 50;
Phire = 46;
Phirn = 6;
Phite = 46.1;
Phitn = 6.1;
Hrg = 10;
Htg = 10;
Grx = 0;
Gtx = 0;
pol = 1;
dct = 500;
dcr = 500;
press = 1013;
temp = 20;

count = 1;
imin = 20;
for ii = imin:length(d)
    dd = d(1:ii);
    hh = h(1:ii);
    zz = z(1:ii);

        L0(count) = tl_p452(GHz, Tpc, dd, hh, zz, Htg, Hrg, Phite, Phitn, Phire, Phirn, Gtx, Grx, pol, dct, dcr, press, temp, 0);
    L1(count) = tl_p452(GHz, Tpc, dd, hh, zz, Htg, Hrg, Phite, Phitn, Phire, Phirn, Gtx, Grx, pol, dct, dcr, press, temp, 1);
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