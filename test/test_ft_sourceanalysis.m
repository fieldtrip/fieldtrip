function [varargout] = test_ft_sourceanalysis(datainfo, writeflag, version, diagnosticsflag)

% MEM 4gb
% WALLTIME 04:30:00

% TEST test_ft_sourceanalysis 
% TEST ft_sourceanalysis ref_datasets

% writeflag determines whether the output should be saved to disk
% version determines the output directory
% diagnosticsflag determines whether some output will be generated for
% diagnostics, rather than exiting on error

if nargin<1
  datainfo = ref_datasets;
end
if nargin<2
  writeflag = 0;
end
if nargin<3
  version = 'latest';
end
if nargin<4
  diagnosticsflag = 0;
end
diagnostics = {};

% make vol
vol = [];
vol.o = [0 0 4];
vol.r = 12;
vol.unit = 'cm';
vol.type = 'singlesphere';

% 3D folded cortical sheet
load(dccnpath('/home/common/matlab/fieldtrip/data/test/corticalsheet.mat'));
sourcemodel_sheet = [];
sourcemodel_sheet.pos = corticalsheet.pnt(1:100,:);          % FIXME reduce the size of the mesh
sourcemodel_sheet.inside = 1:size(sourcemodel_sheet.pos,1);  % FIXME this should not be needed
sourcemodel_sheet.unit = 'cm';

% 3D regular grid
sourcemodel_grid = [];
sourcemodel_grid.resolution = 2.5;
sourcemodel_grid.xgrid = 'auto';
sourcemodel_grid.ygrid = 'auto';
sourcemodel_grid.zgrid = 'auto';

% small number of dipoles, i.e. regions of interest
sourcemodel_roi = [];
sourcemodel_roi.pos = [0 0 5; 1 0 5; -1 0 5; 0 1 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following section can be used to generate all combinations 
% for the test computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if false

set1 = {
  'sheet'
  'grid'
  'roi'
};

set2 = {
  'freq_mtmfft_fourier_trl'
  'freq_mtmconvol_fourier_trl'
  'freq_mtmfft_trl'
  'freq_mtmfft'
  'freq_mtmconvol_trl'
  'freq_mtmconvol'
  'timelock'
  'timelock_trl'
  'timelock_cov'
  'timelock_cov_trl'
};

set3 = {
  'DICS_keepall'
  'DICS_keepall_rawtrial'
  'DICS_keepnothing'
  'DICS_keepnothing_rawtrial'
  'DICS_refdip'
  'DICS_refchan'
  'DICS_realfilter'
  'DICS_fixedori'
  'MNE_keepall'
  'MNE_keepnothing'
  'MNE_keepall_rawtrial'
  'MNE_keepnothing_rawtrial'
  'LCMV_keepall'
  'LCMV_keepnothing'
  'LCMV_keepall_rawtrial'
  'LCMV_keepnothing_rawtrial'
  'DICS_keepall'
  'DICS_keepall_rawtrial'
  'DICS_keepnothing'
  'DICS_keepnothing_rawtrial'
  'DICS_refdip'
  'DICS_refchan'
  'DICS_realfilter'
  'DICS_fixedori'
  'PCC_keepall'
  'PCC_keepall_rawtrial'
  'PCC_keepnothing'
  'PCC_keepnothing_rawtrial'
  'PCC_refdip'
};

i = 1;
n = length(set1)*length(set2)*length(set3)
for i1=1:length(set1)
for i2=1:length(set2)
for i3=1:length(set3)
combination{i,1} = set1{i1};
combination{i,2} = set2{i2};
combination{i,3} = set3{i3};
i = i + 1;
end
end
end

% it requires manual intervention to update the list below
keyboard

end % constructing all combinations

% this script (in test/private) generates the list of allowed combinations
test_ft_sourceanalysis_combinations_allowed;

% TODO: there's an equivalent list of forbidden combinations that should be
% tested to generate an explicit error -> using freq data for time domain
% methods, or using tlck data for freq domain methods.

type = {datainfo.type}'
sel  = strcmp(type, 'meg');
datainfo = datainfo(sel);

for k = 1:numel(datainfo)
    
    %%%%%%FIXME added by JM for the time being, neuromag306 with both a
    %%%%%%grad and an elec fails
    if strcmp(datainfo(k).datatype, 'neuromag306')
        removeelec = true;
    else
        removeelec = false;
    end
    
for j = 1:size(combination,1)

  clear timelock freq data

  sourcemodel        = combination{j,1};
  datarepresentation = combination{j,2};
  algorithm          = combination{j,3};

  switch sourcemodel
  case 'sheet'
    grid = sourcemodel_sheet;
  case 'grid'
    grid = sourcemodel_grid;
  case 'roi'
    grid = sourcemodel_roi;
  end

  switch datarepresentation(1)   % this starts with timelock or freq
  case 'f'
    inputfile = fullfile(datainfo(k).origdir,version,'freq',    datainfo(k).type,[datarepresentation '_' datainfo(k).datatype '.mat']);
    load(inputfile);
    data = freq;
    sourcerepresentation = ['source_' sourcemodel '_' datarepresentation(6:end)]; % drop the 'freq' from the name
  case 't'
    inputfile = fullfile(datainfo(k).origdir,version,'timelock',datainfo(k).type,[datarepresentation '_' datainfo(k).datatype '.mat']);
    load(inputfile);
    data = timelock;
    sourcerepresentation = ['source_' sourcemodel '_' datarepresentation];
  end

  if removeelec && isfield(data, 'elec') 
      data = rmfield(data, 'elec');
  end
  testfunction = str2func(sprintf('sourceanalysis_%s', algorithm));

  outputfile = fullfile(datainfo(k).origdir,version,'source',datainfo(k).type,[sourcerepresentation '_' algorithm '_' datainfo(k).datatype '.mat']);
  try
    fprintf('----------------------------------------------------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------------------------------------------------\n');
    fprintf('inputfile  = %s\n', inputfile);
    fprintf('outputfile = %s\n', outputfile);
    fprintf('----------------------------------------------------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------------------------------------------------\n');

    % execute the actual function that performs the computation
      %try
        source = testfunction(data, grid, vol);
      %catch me
      %  if strcmp(me.message, 'method ''mne'' is unsupported for source reconstruction in the frequency domain');
      %    % continue, and ensure 'source' to exist
      %    source = [];
      %    load(outputfile); % REMOVE THIS ONCE THE FUNCTION RUNS THROUGH AGAIN
      %  elseif strcmp(me.message, 'rawtrial in combination with pcc has been temporarily disabled');
      %    if exist('outputfile', 'file'),
      %      load(outputfile);
      %    else
      %      source = [];
      %    end
      %  else
      %    keyboard;
      %  end
      %end
      
      if writeflag
          save(outputfile, 'source');
      else
          sourcenew = source;
          clear source
          load(outputfile); % this contains the previous "source"
          if isfield(sourcenew, 'cfg'), sourcenew = rmfield(sourcenew, 'cfg'); end% these are different, a.o. due to the callinfo
          source    = rmfield(source, 'cfg');
          if ~diagnosticsflag,
              assert(isequalwithequalnans(source, sourcenew), sprintf('assertion failed: the computed data are different from the data in file %s',outputfile));
          else
              diagnostics{k,j,1} = combination(j,:);
              if isfield(source, 'avg')
                  if issubfield(source, 'avg.pow')
                      diagnostics{k,j,2} = 'pow';
                      diagnostics{k,j,3} = max(source.avg.pow-sourcenew.avg.pow)./max(source.avg.pow);
                  elseif isfield(source.avg, 'mom')
                      diagnostics{k,j,2} = 'mom';
                      diagnostics{k,j,3} = max(source.avg.mom{1}-sourcenew.avg.mom{1})./max(source.avg.mom{1});
                  elseif isfield(source.avg, 'csd')
                      diagnostics{k,j,2} = 'csd';
                      diagnostics{k,j,3} = max(source.avg.csd{1}-sourcenew.avg.csd{1})./max(source.avg.csd{1});
                  end
              elseif isfield(source, 'trial')
                  if isfield(source, 'trial.pow')
                      diagnostics{k,j,2} = 'pow';
                      diagnostics{k,j,3} = max(source.trial(1).pow-sourcenew.trial(1).pow)./max(source.trial(1).pow);
                  elseif isfield(source.avg, 'mom')
                      diagnostics{k,j,2} = 'mom';
                      diagnostics{k,j,3} = max(source.trial(1).mom{1}-sourcenew.trial(1).mom{1})./max(source.trial(1).mom{1});
                  elseif isfield(source.trial(1), 'csd')
                      diagnostics{k,j,2} = 'csd';
                      diagnostics{k,j,3} = max(source.trial(1).csd{1}-sourcenew.trial(1).csd{1})./max(source.trial(1).csd{1});
                  end
              end
          end
      end

  catch me
      
    if strcmp(me.message, sprintf('assertion failed: the computed data are different from the data in file %s',outputfile))
      error('assertion failed: the computed data are different from the data in file %s',outputfile);
    elseif strcmp(me.message, 'method ''mne'' is unsupported for source reconstruction in the frequency domain')
      warning(me.message);
    elseif strcmp(me.message, 'rawtrial in combination with pcc has been temporarily disabled')
      warning(me.message);
    else
      % not all combinations are going to work, give a warning if it fails
      %warning('failed on %s', outputfile);
      error(me.message);
    end
  end

end % combination
end % datainfo

if nargout
  varargout{1} = diagnostics;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MNE subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source = sourceanalysis_MNE_keepall(data, grid, vol)
cfg                   = [];
cfg.channel           = 'MEG';
cfg.method            = 'mne';
cfg.mne.keepleadfield = 'yes';
cfg.mne.keepfilter    = 'yes';
cfg.mne.lambda        = 1e4;
cfg.vol               = vol;
cfg.grid              = grid;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_MNE_keepnothing(data, grid, vol)
cfg                   = [];
cfg.channel           = 'MEG';
cfg.method            = 'mne';
cfg.mne.keepleadfield = 'no';
cfg.mne.keepfilter    = 'no';
cfg.mne.lambda        = 1e4;
cfg.vol               = vol;
cfg.grid              = grid;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_MNE_keepall_rawtrial(data, grid, vol) 
% construct the average spatial filter
source = sourceanalysis_MNE_keepall(data, grid, vol);
% project all trials through the average spatial filter
cfg                   = [];
cfg.channel           = 'MEG';
cfg.method            = 'mne';
cfg.mne.lambda        = 1e4;
cfg.mne.keepleadfield = 'yes';
cfg.mne.keepfilter    = 'yes';
cfg.vol               = vol;
cfg.grid              = grid;
cfg.rawtrial          = 'yes';
cfg.grid.filter       = source.avg.filter;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_MNE_keepnothing_rawtrial(data, grid, vol)
% construct the average spatial filter
source = sourceanalysis_MNE_keepall(data, grid, vol);
% project all trials through the average spatial filter
cfg                   = [];
cfg.method            = 'mne';
cfg.channel           = 'MEG';
cfg.mne.lambda        = 1e4;
cfg.mne.keepleadfield = 'no';
cfg.mne.keepfilter    = 'no';
cfg.vol               = vol;
cfg.grid              = grid;
cfg.rawtrial          = 'yes';
cfg.grid.filter       = source.avg.filter;
source = ft_sourceanalysis(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LCMV subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source = sourceanalysis_LCMV_keepall(data, grid, vol)
cfg                    = [];
cfg.channel            = 'MEG';
cfg.method             = 'lcmv';
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.keepcov       = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.vol                = vol;
cfg.grid               = grid;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_LCMV_keepnothing(data, grid, vol)
cfg                    = [];
cfg.channel            = 'MEG';
cfg.method             = 'lcmv';
cfg.lcmv.keepleadfield = 'no';
cfg.lcmv.keepfilter    = 'no';
cfg.lcmv.keepcov       = 'no';
cfg.lcmv.lambda        = '5%';
cfg.vol                = vol;
cfg.grid               = grid;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_LCMV_keepall_rawtrial(data, grid, vol)
% construct the average spatial filter
source = sourceanalysis_LCMV_keepall(data, grid, vol);
% project all trials through the average spatial filter
cfg                    = [];
cfg.channel            = 'MEG';
cfg.method             = 'lcmv';
cfg.lcmv.keepleadfield = 'yes';
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.keepcov       = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.vol                = vol;
cfg.grid               = grid;
cfg.rawtrial           = 'yes';
cfg.grid.filter        = source.avg.filter;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_LCMV_keepnothing_rawtrial(data, grid, vol)
% construct the average spatial filter
source = sourceanalysis_LCMV_keepall(data, grid, vol);
% project all trials through the average spatial filter
cfg                    = [];
cfg.channel            = 'MEG';
cfg.method             = 'lcmv';
cfg.lcmv.keepleadfield = 'no';
cfg.lcmv.keepfilter    = 'no';
cfg.lcmv.keepcov       = 'no';
cfg.lcmv.lambda        = '5%';
cfg.vol                = vol;
cfg.grid               = grid;
cfg.rawtrial           = 'yes';
cfg.grid.filter        = source.avg.filter;
source = ft_sourceanalysis(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DICS subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function source = sourceanalysis_DICS_keepall(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.keeptrials   = 'yes';
cfg.dics.projectnoise = 'yes';
cfg.dics.feedback = 'none';
% dics options
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.keepmom       = 'yes';
cfg.dics.lambda        = '5%';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_keepall_rawtrial(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.keeptrials   = 'yes';
cfg.dics.projectnoise = 'yes';
cfg.dics.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.keepmom       = 'yes';
cfg.dics.lambda        = '5%';
% get filter
tmp = sourceanalysis_DICS_keepall(data, grid, vol);
cfg.rawtrial    = 'yes';
cfg.grid.filter = tmp.avg.filter;
cfg.feedback    = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_keepnothing(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.keeptrials   = 'no';
cfg.dics.projectnoise = 'no';
cfg.dics.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'no';
cfg.dics.keepleadfield = 'no';
cfg.dics.keepcsd       = 'no';
cfg.dics.lambda        = '5%';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_keepnothing_rawtrial(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.keeptrials   = 'no';
cfg.dics.projectnoise = 'no';
cfg.dics.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'no';
cfg.dics.keepleadfield = 'no';
cfg.dics.keepcsd       = 'no';
cfg.dics.lambda        = '5%';
% get filter
tmp = sourceanalysis_DICS_keepall(data, grid, vol);
cfg.rawtrial    = 'yes';
cfg.grid.filter = tmp.avg.filter;
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_refdip(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.dics.refdip       = [2 5 9];
cfg.dics.keepcsd      = 'yes';
cfg.dics.feedback     = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_refchan(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
if numel(data.label)<40
    cfg.refchan = data.label{5};
else
    cfg.refchan = data.label{40};
end
cfg.dics.keepcsd      = 'yes';
cfg.dics.feedback     = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_realfilter(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.dics.realfilter   = 'yes';
cfg.dics.keepfilter   = 'yes';
cfg.dics.feedback     = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_DICS_fixedori(data, grid, vol)
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.latency      = 0.5;
cfg.dics.fixedori     = 'yes';
cfg.dics.feedback     = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_PCC_keepall(data, grid, vol)
% do PCC
cfg = [];
cfg.method            = 'pcc';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.keepleadfield = 'yes';
cfg.pcc.keepcsd       = 'yes';
cfg.pcc.keepmom       = 'yes';
cfg.pcc.lambda        = '5%';
cfg.channel           = 'MEG';
cfg.frequency         = 10;
cfg.latency           = 0.5;
cfg.vol               = vol;
cfg.grid              = grid;
source                = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_PCC_keepall_rawtrial(data, grid, vol)
% do PCC
cfg = [];
cfg.method            = 'pcc';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.keepleadfield = 'yes';
cfg.pcc.keepcsd       = 'yes';
cfg.pcc.keepmom       = 'yes';
cfg.pcc.lambda        = '5%';
cfg.channel           = 'MEG';
cfg.frequency         = 10;
cfg.latency           = 0.5;
cfg.vol               = vol;
cfg.grid              = grid;
tmp                   = ft_sourceanalysis(cfg, data);
% get filter
tmp = sourceanalysis_PCC_keepall(data, grid, vol);
cfg.rawtrial    = 'yes';
cfg.grid.filter = tmp.avg.filter;
cfg.pcc.feedback    = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_PCC_keepnothing(data, grid, vol)
% do PCC
cfg = [];
cfg.method            = 'pcc';
cfg.pcc.keepfilter    = 'no';
cfg.pcc.keepleadfield = 'no';
cfg.pcc.keepcsd       = 'no';
cfg.pcc.keepmom       = 'no';
cfg.pcc.lambda        = '5%';
cfg.channel           = 'MEG';
cfg.frequency         = 10;
cfg.latency           = 0.5;
cfg.vol               = vol;
cfg.grid              = grid;
source                = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_PCC_keepnothing_rawtrial(data, grid, vol)
% do PCC
cfg = [];
cfg.method            = 'pcc';
cfg.pcc.keepfilter    = 'no';
cfg.pcc.keepleadfield = 'no';
cfg.pcc.keepcsd       = 'no';
cfg.pcc.keepmom       = 'no';
cfg.pcc.lambda        = '5%';
cfg.channel           = 'MEG';
cfg.frequency         = 10;
cfg.latency      = 0.5;
cfg.vol               = vol;
cfg.grid              = grid;
tmp                   = ft_sourceanalysis(cfg, data);
% get filter
tmp = sourceanalysis_PCC_keepall(data, grid, vol);
cfg.rawtrial    = 'yes';
cfg.grid.filter = tmp.avg.filter;
cfg.pcc.feedback    = 'none';
source = ft_sourceanalysis(cfg, data);

function source = sourceanalysis_PCC_refdip(data, grid, vol)
% do PCC
cfg = [];
cfg.method            = 'pcc';
cfg.channel           = 'MEG';
cfg.frequency         = 10;
cfg.latency      = 0.5;
cfg.vol               = vol;
cfg.grid              = grid;
cfg.refdip            = [2 5 9];
%cfg.keepcsd           = 'yes';   % keepcsd is ALWAYS ON with PCC
source                = ft_sourceanalysis(cfg, data);
