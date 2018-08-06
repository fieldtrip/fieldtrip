function test_bug2222

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqcomparison ft_math

% the bug reported by Izabela triggered a short discussion, concluding that
% ft_freqcomparison should be deprecated. This test script tests the (old)
% function and ensures that ft_math (it's replacement) performs the same or better.

freq1 = [];
freq1.dimord = 'chan_freq';
freq1.freq = 1:10;
freq1.label = {'1', '2', '3'};
freq1.powspctrm = randn(3,10);
freq1.powspctrm(1) = 1;

freq2 = [];
freq2.dimord = 'chan_freq';
freq2.freq = 1:10;
freq2.label = {'1', '2', '3'};
freq2.powspctrm = randn(3,10);
freq2.powspctrm(1) = 2;

cfg = [];
cfg.comparisontype = 'absolute';
freq3a = ft_freqcomparison(cfg, freq1, freq2); % this is 2 - 1
cfg.comparisontype = 'relative';
freq3b = ft_freqcomparison(cfg, freq1, freq2); % this is 2 / 1
cfg.comparisontype = 'relchange';
freq3c = ft_freqcomparison(cfg, freq1, freq2); % this is (2-1)/1

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'subtract';
freq4a = ft_math(cfg, freq2, freq1);  % this is 2 - 1
cfg.operation = 'divide';
freq4b = ft_math(cfg, freq2, freq1);  % this is 2 / 1
cfg.operation = 'divide';
freq4c = ft_math(cfg, freq4a, freq1); % this is (2-1)/1 in two steps

% it turned out that the output of ft_freqcomparison did not contain a cfg in this
% case, but also that it would contain an incorrect cfg if the input has a cfg. I
% fixed that along with making this test script.

assert(isequal(rmfield(freq3a, 'cfg'), rmfield(freq4a, 'cfg')), 'failed for a');
assert(isequal(rmfield(freq3b, 'cfg'), rmfield(freq4b, 'cfg')), 'failed for b');
assert(isequal(rmfield(freq3c, 'cfg'), rmfield(freq4c, 'cfg')), 'failed for c');
