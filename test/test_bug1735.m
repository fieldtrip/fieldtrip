function test_bug1735

% TEST test_bug1735
% TEST ft_appenddata

load test_bug1735.mat

cfg = [];
test1 = ft_appenddata([], data_trgtOnstimOnBTrgtLDetected_sourceL, data_trgtOnstimOnBTrgtLNotDetected_sourceL);

assert(length(test1.trial)==79);
