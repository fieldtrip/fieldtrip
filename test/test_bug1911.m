function test_bug1911

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_databrowser
% DATA private

% When ft_movieplotER is called within ft_databrowser by right-clicking on a
% segment of data and selecting ft_movieplotER the movieplot is opened in the
% ft_databrowser window. Only a quarter of the topoplot is visible in this case. 
%
% I have only experienced this problem in Windows. 

cd(dccnpath('/project/3031000.02/test'));
load bug1911.mat

cfg = [];
cfg.layout = 'easycapM10.lay';
ft_databrowser(cfg,timelock);

