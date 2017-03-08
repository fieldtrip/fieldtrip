function test_bug1735

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_appenddata

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1735.mat

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

% following up with the databrowser, it turned out that there was a problem
% it is due to combined.sampleinfo having overlapping trials

try
  ft_fetch_data(combined, 'begsample', 1, 'endsample', 50, 'allowoverlap', false);
  status = false; % an error should be given by ft_fetch_data
catch
  status = true;
end
assert(status, 'ft_fetch_data did not detect the error properly');






