function test_bug3397

% WALLTIME 00:10:00
% MEM 1gb
% TEST ft_appenddata ft_preamble_init

data       = [];
data.label = {'chan01';'chan02';'chan03'};
for k = 1:10
trial{k} = randn(3,100);
time{k}  = (0:99)./100;
end
data.time  = time;
data.trial = trial;
 
alldata{1} = data;
alldata{2} = data;
alldata{3} = data;

tmp = ft_appenddata([], alldata{:});
assert(numel(tmp.trial)==30);

try
  ft_appenddata(alldata{:});
catch ME
  if strncmp(ME.message,'It seems as if the first', 24)
    % this is the expected error
  else
    rethrow(ME);
  end
end


