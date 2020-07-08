function test_ft_timelockanalysis(datainfo, writeflag, version)

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

% ft_timelockanalysis ref_datasets

% writeflag determines whether the output should be saved to disk
% version determines the output directory

if nargin<1 || isempty(datainfo)
  datainfo = ref_datasets;
end
if nargin<2 || isempty(writeflag)
  writeflag = 0;
end
if nargin<3
  version = 'latest';
end

for k = 1:numel(datainfo)
  disp(['now you are in k=' num2str(k)]);
  
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'yes', 'yes');
  fname   = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,['timelock_cov_trl_',datainfo(k).datatype]);
  comparedata(datanew, fname);
  
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'yes', 'no');
  fname   = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,['timelock_cov_',datainfo(k).datatype]);
  comparedata(datanew, fname);
  
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'no', 'yes');
  fname   = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,['timelock_trl_',datainfo(k).datatype]);
  comparedata(datanew, fname);
  
  datanew = timelockanalysis10trials(datainfo(k), writeflag, version, 'no', 'no'); % should be the latest
  fname   = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,['timelock_',datainfo(k).datatype]);
  comparedata(datanew, fname);
  
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

function comparedata(datanew, fname)

tmp = load(fname);
if isfield(tmp, 'data')
  data = tmp.data;
elseif isfield(tmp, 'datanew')
  data = tmp.datanew;
elseif isfield(tmp, 'timelock')
  data = tmp.timelock;
end

try
  % these are sensitive to numerical erorrs, especially for MEG data and when the covariance is close to zero
  datanew.cov = trimnumericalerror(datanew.cov, eps*1e3);
  data.cov    = trimnumericalerror(data.cov   , eps*1e3);
end

datanew = removefields(datanew, 'cfg'); % these are per construction different if writeflag = 0
data    = removefields(data,    'cfg');
[ok,msg] = isalmostequal(data, datanew, 'reltol', eps*1e8);
if ~ok
  error('stored and computed data not identical: %s', msg{:});
end
