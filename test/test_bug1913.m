function test_bug1913

% MEM 2gb

% TEST test_bug1913

fn1 = '/home/common/matlab/fieldtrip/data/test/original/eeg/eeprobe/martlbc1.cnt';
fn2 = '/home/common/matlab/fieldtrip/data/test/original/eeg/eeprobe/martlbc1s1.avr';
fn3 = '/home/common/matlab/fieldtrip/data/test/original/eeg/eeprobe/VP_05_1.cnt';

hdr1 = ft_read_header(fn1);
hdr2 = ft_read_header(fn2);
hdr3 = ft_read_header(fn3);

dat1 = ft_read_data(fn1); clear dat*
dat2 = ft_read_data(fn2); clear dat*
dat3 = ft_read_data(fn3); clear dat*
