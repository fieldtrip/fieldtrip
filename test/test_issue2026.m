function test_issue2026

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_read_event

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/issue2026/yev0004_New-Coherent.vmrk')

%%

% the issue is that the code that reads the markers originally uses strsplit for tokenization. The code expects empty spaces between the splitting tokens to be empty elements in the output of the splitting, but strsplit does not do this. As opposed to FT's tokenize function, or the Mathworks recommend split function (the suggested alternative split() is only present in MATLAB 2016b onwards, so strictly speaking it could be used.

event = ft_read_event(filename);
val   = {event.value}';
if any(cellfun(@isnumeric, val) & ~cellfun(@isempty, val))
  % this means that the sample index ended up unexpectedly at the level of the value, which is incorrect
  error('something went wrong with the string splitting');
end 
