function test_bug3397

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_appenddata ft_preamble_init

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
  failed = false;  
catch
  failed = true;
end
assert(failed, 'the expected error was not thrown');


