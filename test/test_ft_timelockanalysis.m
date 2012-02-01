function test_ft_timelockanalysis(datainfo, writeflag)

% TEST test_ft_timelockanalysis
% ft_timelockanalysis test_datasets

% the optional writeflag determines whether the output should be saved
% to disk

if nargin<2
  writeflag = 0;
end
if nargin<1
  datainfo = test_datasets;
end
for k = 1:numel(datainfo)
  datanew = timelockanalysis10trials(datainfo(k), writeflag);
  
  fname = [datainfo(k).origdir,'timelock/',datainfo(k).type,'timelock_',datainfo(k).datatype];
  tmp = load(fname);
  if isfield(tmp, 'data')
    data = tmp.data;
  elseif isfield(tmp, 'datanew')
    data = tmp.datanew;
  else isfield(tmp, 'timelock')
    data = tmp.timelock;
  end
  
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  data    = rmfield(data,    'cfg');
  assert(isequalwithequalnans(data, datanew));
end

function [tlck] = timelockanalysis10trials(dataset, writeflag)

cfg        = [];
cfg.inputfile  = [dataset.origdir,'raw/',dataset.type,'preproc_',dataset.datatype];
if writeflag
  cfg.outputfile = [dataset.origdir,'timelock/',dataset.type,'timelock_',dataset.datatype];
end
tlck = ft_timelockanalysis(cfg);
