function test_ft_timelockanalysis(datainfo, writeflag, version)

% MEM 1500mb
% WALLTIME 00:10:00

% ft_timelockanalysis ref_datasets

% writeflag determines whether the output should be saved to disk
% version determines the output directory

if nargin<1
  datainfo = ref_datasets;
end
if nargin<2
  writeflag = 0;
end
if nargin<3
  version = 'latest';
end

for k = 1:numel(datainfo)
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'yes', 'yes');
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'yes', 'no');
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'no', 'yes');
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'no', 'no'); % should be the latest
  
  fname = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,['timelock_',datainfo(k).datatype]);
  tmp = load(fname);
  if isfield(tmp, 'data')
    data = tmp.data;
  elseif isfield(tmp, 'datanew')
    data = tmp.datanew;
  else isfield(tmp, 'timelock')
    data = tmp.timelock;
  end
  
  datanew = removefields(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  data    = removefields(data,    'cfg');
  [ok,msg] = isalmostequal(data, datanew,'reltol',eps*1e6);
  disp(['now you are in k=' num2str(k)]);
  if ~ok
    error('stored and computed data not identical: %s', msg{:});
  end
end

function [timelock] = timelockanalysis10trials(dataset, writeflag, version, covariance, keeptrials)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_timelockanalysis') && exist('timelockanalysis')
  eval('ft_timelockanalysis = @timelockanalysis;');
end

if isempty(covariance)
  covariance = 'no';
  % covariance = 'yes';
end

if isempty(keeptrials)
  keeptrials = 'no';
  % keeptrials = 'yes';
end

% the file names should distinguish between the cfg.covariance and cfg.keeptrials option
postfix = '';
switch covariance
case 'no'
  % don't change
case 'yes'
  postfix = [postfix 'cov_'];
otherwise
  error('unexpected keeptrials');
end

% the file names should distinguish between the cfg.covariance and cfg.keeptrials option
switch keeptrials
case 'no'
  % don't change
case 'yes'
  postfix = [postfix 'trl_'];
otherwise
  error('unexpected keeptrials');
end

cfg = [];
cfg.keeptrials = keeptrials;
cfg.covariance = covariance;
cfg.inputfile  = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype]);
if writeflag
  cfg.outputfile = fullfile(dataset.origdir,version,'timelock',dataset.type,['timelock_',postfix,dataset.datatype]);
end

if ~strcmp(version, 'latest') && str2num(version)<20100000
  % -- HISTORICAL --- older FieldTrip versions don't support inputfile and outputfile
  load(cfg.inputfile, 'data');
  timelock = ft_timelockanalysis(cfg, data);
  save(cfg.outputfile, 'timelock');
else
  timelock = ft_timelockanalysis(cfg);
end

