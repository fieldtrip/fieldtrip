function test_ft_componentanalysis(datainfo, writeflag, version)

% MEM 2gb
% WALLTIME 00:10:00

% ft_componentanalysis ref_datasets

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
  datanew = componentanalysis(datainfo(k), writeflag, version);
  
  fname = fullfile(datainfo(k).origdir,version,'comp',datainfo(k).type,['comp_',datainfo(k).datatype]);
  tmp = load(fname);
  if isfield(tmp, 'comp')
    data = tmp.comp;
  end
  
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  data    = rmfield(data,    'cfg');
  
  % if data is rank-deficient, the last columns of the mixing/unmixing
  % matrices are arbitrary, and should thus not be compared
  rankDiff = size(data.trial{1},1) - rank(data.trial{1});
  if rankDiff == size(data.trial{1},1)
    % massive rank deficiency (i.e., identical data in all channels)
    % best to just not test the mixing matrices at all, just surrogate test
    % data
    data = rmfield(data, 'topo');
    data = rmfield(data, 'unmixing');
    datanew = rmfield(datanew, 'topo');
    datanew = rmfield(datanew, 'unmixing');
  elseif rankDiff > 0
    data.topo(:,end-rankDiff:end) = 0;
    data.unmixing(end-rankDiff:end,:) = 0;
    datanew.topo(:,end-rankDiff:end) = 0;
    datanew.unmixing(end-rankDiff:end,:) = 0;
  end
  
  [ok, msg] = isalmostequal(data, datanew, 'abstol', 1e-5, 'diffabs', 1);
  disp(['now you are in k=' num2str(k)]);
  if ~ok
    disp(msg);
    error('there were differences between reference and new data, see above for details');
  end
end

function [comp] = componentanalysis(dataset, writeflag, version)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_componentanalysis') && exist('componentanalysis')
  eval('ft_componentanalysis = @componentanalysis;');
end

cfg = [];
cfg.method = 'pca';

switch dataset.datatype
  case {'bti148' 'bti248' 'bti248grad' 'ctf151' 'ctf275' 'ctf64' 'itab153'  'neuromag122' 'neuromag306' 'yokogawa160'}
    cfg.channel = 'MEG';
  otherwise
    cfg.channel = 'all';
end

cfg.inputfile  = fullfile(dataset.origdir,version,'raw',dataset.type,['preproc_',dataset.datatype]);
outputfile     = fullfile(dataset.origdir,version,'comp',dataset.type,['comp_',dataset.datatype]);
if writeflag
  cfg.outputfile = outputfile;
end

if ~strcmp(version, 'latest') && str2double(version)<20100000
  % -- HISTORICAL --- older FieldTrip versions don't support inputfile and outputfile
  try
    % use the previous random seed
    load(outputfile, 'comp');
    if isfield(comp.cfg.callinfo, 'randomseed')
      cfg.randomseed = comp.cfg.callinfo.randomseed;
    end
  catch
    % use a new random seed
  end
  
  load(cfg.inputfile, 'data');
  comp = ft_componentanalysis(cfg, data);
  save(cfg.outputfile, 'comp');
else
  try
    % use the previous random seed
    load(outputfile, 'comp');
    if isfield(comp.cfg.callinfo, 'randomseed')
      cfg.randomseed = comp.cfg.callinfo.randomseed;
    end
  catch
    % use a new random seed
  end
  
  comp = ft_componentanalysis(cfg);
end

