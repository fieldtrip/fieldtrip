function test_bug1924

% MEM 1500mb
% WALLTIME 00:10:00

% TEST read_deymed_dat read_deymed_ini ft_filetype ft_read_header ft_read_data ft_read_event

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1924'));

datfile = 'R2.09.Dat';
edffile = 'R2.09-export.edf';

hdr1 = ft_read_header(datfile);
hdr2 = ft_read_header(edffile);

dat1 = ft_read_data(datfile);
dat2 = ft_read_data(edffile);

sel1 = (1:25600);
sel2 = 1:25600;

x = dat1(2,sel1);
y = dat2(2,sel2);

figure
plot(x, 'b')
hold on
plot(y, 'r.')

z = xcorr(x, y, 'coeff');
[m, i] = max(z);

% this should peak in the middle
figure
plot(z)

assert(isalmostequal(x, y, 'reltol', 0.1));
% They are quite different, hence the tolerance of 10%. This might well be
% due to the BESA conversion, in which the numbers have to be rescaled to
% fit in the EDF file. On average the calibration seems quite OK.

