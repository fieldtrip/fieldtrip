function test_bug2163

% TEST test_bug2163
% TEST ft_read_spike read_neuralynx_nse

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2163/test.nse');
spike = ft_read_spike(filename);
assert(size(spike.waveform{1},2)==64); % there should be 64 samples per spike snippet in this test file

filename = dccnfilename('/home/common/matlab/fieldtrip/data/test/bug2163/Sc1.nse');
spike = ft_read_spike(filename);
assert(size(spike.waveform{1},2)==32); % this is an original file with 32 samples