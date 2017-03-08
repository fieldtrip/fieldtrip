function test_bug1878

% MEM 3gb
% WALLTIME 00:10:00

% TEST ft_artifact_clip

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'));
load bug1878.mat

% this is how I was able to reproduce it

% restoredefaultpath
% addpath(dccnpath('/home/common/matlab/fieldtrip-20120630'));
% [cfg1, artifact1] = ft_artifact_clip(cfg, data1);
% 
% restoredefaultpath
% addpath(dccnpath('/home/common/matlab/fieldtrip'));
% [cfg2, artifact2] = ft_artifact_clip(cfg, data1);

% from now on the error should not happen any more with any future version of fieldtrip
[cfg2, artifact2] = ft_artifact_clip(cfg, data1);

assert(isempty(artifact2), 'this dataset contains no clipping artifacts');

