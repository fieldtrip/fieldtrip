function test_bug1913

% WALLTIME 00:10:00
% MEM 1500mb


fn1 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/eeprobe/martlbc1.cnt');
fn2 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/eeprobe/martlbc1s1.avr');
fn3 = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/eeprobe/VP_05_1.cnt');

hdr1 = ft_read_header(fn1);
hdr2 = ft_read_header(fn2);
hdr3 = ft_read_header(fn3);

dat1 = ft_read_data(fn1); clear dat*
dat2 = ft_read_data(fn2); clear dat*
dat3 = ft_read_data(fn3); clear dat*
