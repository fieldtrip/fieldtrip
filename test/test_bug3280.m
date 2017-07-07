function test_bug3280

% WALLTIME 00:10:00
% MEM 2gb

%%
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3280'));

dataset = 'A1826_comparisonSEF_20161212_02.ds';

%%

hdr = ft_read_header(dataset);

grad = ft_read_sens(dataset, 'senstype', 'meg');
elec = ft_read_sens(dataset, 'senstype', 'eeg');

assert(isfield(grad,'coilpos'));
assert(isfield(elec,'elecpos'));


