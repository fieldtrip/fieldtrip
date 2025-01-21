function test_bug2163

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_spike read_neuralynx_nse
% DATA private

filename = dccnpath('/project/3031000.02/test/bug2163/test.nse');
spike = ft_read_spike(filename);
assert(size(spike.waveform{1},2)==64); % there should be 64 samples per spike snippet in this test file

filename = dccnpath('/project/3031000.02/test/bug2163/Sc1.nse');
spike = ft_read_spike(filename);
assert(size(spike.waveform{1},2)==32); % this is an original file with 32 samples
