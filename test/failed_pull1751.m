function test_pull1751

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_nirs_signalqualityindex

% load sample set of data
fname = dccnpath('/home/common/matlab/fieldtrip/data/test/original/nirs/artinis/Sofia/sample_data.oxy4');

cfg = [];
cfg.dataset = fname; 
cfg.channel = 'nirs';
data_test = ft_preprocessing(cfg);

% load expected data output
load('/home/common/matlab/fieldtrip/data/test/sqi/output_sqi.mat');
 
% compute sqi
cfg = [];
data_out = ft_nirs_signalqualityindex(cfg, data_test);

% check whether it matches the expected output
assert(isequaln(data_out.trial, output_trials))

% check whether the size of the output matrix is consistent with the input one
assert(isequal(size(output_trials{1}),size(data_test.trial{1}))) 

%%%% check whether it raises error when there are problems with the labels
% case 1: channel labels do not match the expected order
data_non_ordered = data_test;
labels = data_non_ordered.label;
idx_rand = randperm(numel(labels));
labels_new = labels(idx_rand);

trials = data_non_ordered.trial;
trials_new{1, 1} = trials{1, 1}(idx_rand,:);

data_non_ordered.trial = trials_new;
data_non_ordered.label = labels_new;

try
    ft_nirs_signalqualityindex(cfg, data_non_ordered);
    raises_error = 0;
catch 
    raises_error = 1;
end

if ~raises_error
    error('Does not raise error when labels do not match expected order.')
end
    
% case 2: more than one label associated to the same channel
data_wrong_num_labels = data_test;
labels = data_wrong_num_labels.label;
labels{3} = 'Rx1-Tx1 [900nm]';
data_wrong_num_labels.label = labels;

try
    ft_nirs_signalqualityindex(cfg, data_wrong_num_labels);
    raises_error = 0;
catch 
    raises_error = 1;
end

if ~raises_error
    error('Does not raise error when channel combination is not associated with a number of wavelengths different than 2.')
end
