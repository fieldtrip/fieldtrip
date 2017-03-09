function test_bug97

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_selectdata

% this script tests the solution to bug 97:
% Selectdata is working in an unexpected way. If the user explicitly states
% cfg.trials = [], no trials are to be processed. The current functionality is
% that if isempty(selrpt), no selection of trials takes place, and thus all
% trials are processed subsequently. This is the result of keyval being used, but
% the output of keyval being inappropriately processed. 
% Solution: if the user explicitly specifies a 'key', then use the next input
% argument as the corresponding 'value'. In such case [] has a definite meaning
% (rather than not being defined by the user)

datafile   = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds');

cfg          = [];
cfg.datafile = datafile;
cfg.trl      = [[1001:1000:10001]' [2000:1000:11000]' round(randn(10,1)*100)];
cfg.trl(:,4) = [ones(5,1); ones(5,1)*2];
cfg.continuous = 'yes';
data           = ft_preprocessing(cfg);

% now subselect no trials
cfg        = [];
cfg.trials = [];
datax      = ft_preprocessing(cfg, data); % this reproduces the behaviour

% after adjusting selectdata
datay      = ft_preprocessing(cfg, data);

