function test_bug1811

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_analysispipeline

% loading meg + eyelink channels data appended with function_handle variable
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1811.mat'));
 
cfg              =[];
cfg.feedback     ='no';
[script]=ft_analysispipeline(cfg,meeg);

if isempty(script)
   error('script is empty');
end
clear script;

%check ctf275 meg preoprocessed data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'));
script = ft_analysispipeline(cfg,data);

if isempty(script)
   error('script is empty');
end
