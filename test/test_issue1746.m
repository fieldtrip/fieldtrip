function test_issue1746

% WALLTIME 00:10:00
% MEM 4gb

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1746'));

%%

hdr = ft_read_header('FS001.csv');
dat = ft_read_data('FS001.csv');
evt = ft_read_event('FS001.csv');

%%

sel = evt(1).sample + (-100:100);

plot(sel, dat(4,sel), '.');
hold on
plot(evt(1).sample, evt(1).value, 'go');
plot(evt(2).sample, evt(2).value, 'go');

%%

cfg = [];
cfg.dataset = 'FS001.csv';
data = ft_preprocessing(cfg);

cfg = [];
cfg.event = evt;
cfg.channel = 'position*';
cfg.preproc.demean = 'yes';
cfg.blocksize = 30;
ft_databrowser(cfg, data);


