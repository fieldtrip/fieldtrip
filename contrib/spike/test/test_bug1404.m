function test_bug1404

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_bug1404
% TEST ft_read_header ft_read_data ft_read_spike

% this is only available on the DCCN mentat cluster
dataset = '/home/common/matlab/fieldtrip/data/test/original/neurosim/signals';
hdr = ft_read_header(dataset);

% this should also work
dataset = '/home/common/matlab/fieldtrip/data/test/original/neurosim';
hdr = ft_read_header(dataset);

cfg = [];
cfg.dataset = dataset;
data = ft_preprocessing(cfg);

% the first channel is the time axis
assert(all((data.time{1} - data.trial{1}(1,:))<0.01/hdr.Fs))

spike = ft_read_spike(dataset);

ft_checkdata(data, 'datatype', 'raw');
ft_checkdata(spike, 'datatype', 'spike');
