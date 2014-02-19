function test_bug1270

% MEM 1500mb
% WALLTIME 00:10:00

% script to verify and fix an error related to using a mask in
% ft_singleplotER
%
% ??? Undefined function or variable 'xvector'.
% 
% Error in ==> ft_singleplotER at 477
%    maskdatavector = reshape(mean(datmask,1), [1 numel(xvector)]);

clear all

mask = repmat([0 0 0 1 1 1 1 1 0 0],2,1);

timelock1 = [];
timelock1.label = {'1' '2'};
timelock1.time  = 1:10;
timelock1.dimord = 'chan_time';
timelock1.avg = randn(2,10);
timelock1.mask = timelock1.avg.*mask;

cfg = [];
cfg.channel = '1';
cfg.maskparameter = 'mask';
ft_singleplotER(cfg, timelock1);

% changing xvector into xval does the trick
