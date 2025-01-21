function test_bug3280

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY
% DATA private


%%
cd(dccnpath('/project/3031000.02/test/bug3280'));

dataset = 'A1826_comparisonSEF_20161212_02.ds';

%%

hdr = ft_read_header(dataset);

grad = ft_read_sens(dataset, 'senstype', 'meg');
elec = ft_read_sens(dataset, 'senstype', 'eeg');

assert(isfield(grad,'coilpos'));
assert(isfield(elec,'elecpos'));


