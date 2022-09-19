function test_issue623

% MEM 2gb
% WALLTIME 00:20:00
% DEPENCENCY

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/brainvision/Mischa.vhdr');

ft_hastoolbox('compat/obsolete', 1); % ft_qualitycheck has moved to be obsolete

%%

cfg = [];
cfg.dataset = filename;
cfg.saveplot = 'no';
cfg.savemat = 'no';
[info, timelock, freq, summary] = ft_qualitycheck(cfg);
