function test_bug1403

% MEM 1500mb
% WALLTIME 0:03:01

% TEST test_bug1403
% TEST ft_read_header

cd /home/common/matlab/fieldtrip/data/test/bug1403

cfg=[];
cfg.dataset = 'LauraPP1_SEM_MATCH_Average_AUT2.vhdr';
data = ft_preprocessing(cfg);

cfg=[];
cfg.layout = '61chan_MPI.lay';
cfg.interactive = 'yes';
figure
ft_multiplotER(cfg, data);

figure
plot(data.time{1}, data.trial{1})
legend(data.label)

