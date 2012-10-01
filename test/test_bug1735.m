function test_bug1735

% TEST test_bug1735
% TEST ft_appenddata

load test_bug1735.mat

part1 = data_trgtOnstimOnBTrgtLDetected_sourceL;
part2 = data_trgtOnstimOnBTrgtLNotDetected_sourceL;

cfg = [];
combined = ft_appenddata(cfg, part1, part2);

assert(length(combined.trial)==79);

N1 = size(part1.sampleinfo,1);
N2 = size(part2.sampleinfo,1);

for i=1:N1
  disp(i)
  assert(isequal(part1.trial{i}, combined.trial{i}));
end

for i=1:N2
  disp(i)
  assert(isequal(part2.trial{i}, combined.trial{i+N1}));
end

% the data in the two conditions should be different
for i=1:N1
  for j=1:N2
    assert(~isequal(part1.trial{i}, part2.trial{j}));
  end
end
