function test_bug1811

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_bug1811
% TEST ft_analysispipeline

% loading meg + eyelink channels data appended with function_handle variable
cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1811.mat
 
cfg              =[];
cfg.feedback     ='no';
[script, details]=ft_analysispipeline(cfg,meeg);

if isempty(script) || isempty(details);
   error('script or details output are empty');
end
clear script details;

%check ctf275 meg preoprocessed data
load '/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf275.mat'
[script, details]=ft_analysispipeline(cfg,data);

if isempty(script) || isempty(details);
   error('script or details output are empty');
end
