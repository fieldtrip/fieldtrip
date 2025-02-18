function test_bug1913

% MEM 6gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header ft_read_data

fn1 = dccnpath('/project/3031000.02/test/original/eeg/eeprobe/martlbc1.cnt');
fn2 = dccnpath('/project/3031000.02/test/original/eeg/eeprobe/martlbc1s1.avr');
fn3 = dccnpath('/project/3031000.02/test/original/eeg/eeprobe/VP_05_1.cnt');

hdr1 = ft_read_header(fn1);
hdr2 = ft_read_header(fn2);
hdr3 = ft_read_header(fn3);

dat1 = ft_read_data(fn1); clear dat*
dat2 = ft_read_data(fn2); clear dat*
dat3 = ft_read_data(fn3); clear dat*
