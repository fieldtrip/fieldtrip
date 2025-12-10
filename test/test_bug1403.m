function test_bug1403

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header
% DATA private

cfg=[];
cfg.dataset = dccnpath('/project/3031000.02/test/bug1403/LauraPP1_SEM_MATCH_Average_AUT2.vhdr');
data = ft_preprocessing(cfg);

cfg=[];
cfg.layout = dccnpath('/project/3031000.02/test/bug1403/61chan_MPI.lay');
cfg.interactive = 'yes';
figure
ft_multiplotER(cfg, data);

figure
plot(data.time{1}, data.trial{1})
legend(data.label)

