function test_issue2082

% MEM 1gb
% WALLTIME 00:30:00
% DEPENDENCY ft_checkconfig
% DATA public


cfg = [];
cfg.dataset = string(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds'));
cfg = ft_checkconfig(cfg, 'dataset2files', 'yes');

filename = string(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.shape'));
hs = ft_read_headshape(filename);
