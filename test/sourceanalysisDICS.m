function sourceanalysisDICS(data, grid, vol)

% general source
% sourceanalysisDICS_keepall(data, grid, vol);
% sourceanalysisDICS_keepall_rawtrial(data, grid, vol);
% sourceanalysisDICS_keepnothing(data, grid, vol);
% sourceanalysisDICS_keepnothing_rawtrial(data, grid, vol);

% dics specific
source = sourceanalysisDICS_refdip(data, grid, vol);
sourceanalysisDICS_refchan(data, grid, vol);
sourceanalysisDICS_realfilter(data, grid, vol);
sourceanalysisDICS_fixedori(data, grid, vol);

end

function source = sourceanalysisDICS_keepall(data, grid, vol)
% dataset should be:
% if ispc
%   dataset =
%   H:\common\matlab\fieldtrip\data\test\latest\freq\meg\freq_mtmfft_fourier_trl_ctf275.mat;
% elseif isunix
%   dataset = fullfile('/home', 'common', 'matlab', 'fieldtrip', 'data', 'test', 'latest', 'freq', 'meg', 'freq_mtmfft_fourier_trl_ctf275.mat');
% end

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.keeptrials   = 'yes';
cfg.projectnoise = 'yes';
cfg.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.keepmom       = 'yes';
cfg.dics.lambda        = '5%';

source = ft_sourceanalysis(cfg, data);

end


function source = sourceanalysisDICS_keepall_rawtrial(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.keeptrials   = 'yes';
cfg.projectnoise = 'yes';
cfg.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'yes';
cfg.dics.keepleadfield = 'yes';
cfg.dics.keepcsd       = 'yes';
cfg.dics.keepmom       = 'yes';
cfg.dics.lambda        = '5%';

% get filter
tmp = sourceanalysisDICS_keepall(data, grid, vol);
cfg.rawtrial    = 'yes';
cfg.grid.filter = tmp.avg.filter;
cfg.feedback    = 'none';

source = ft_sourceanalysis(cfg, data);

end

function source = sourceanalysisDICS_keepnothing(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.keeptrials   = 'no';
cfg.projectnoise = 'no';
cfg.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'no';
cfg.dics.keepleadfield = 'no';
cfg.dics.keepcsd       = 'no';
cfg.dics.lambda        = '5%';

source = ft_sourceanalysis(cfg, data);

end

function source = sourceanalysisDICS_keepnothing_rawtrial(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.keeptrials   = 'no';
cfg.projectnoise = 'no';
cfg.feedback     = 'none';
% dics options
cfg.dics.keepfilter    = 'no';
cfg.dics.keepleadfield = 'no';
cfg.dics.keepcsd       = 'no';
cfg.dics.lambda        = '5%';

% get filter
tmp = sourceanalysisDICS_keepall(data, grid, vol);
cfg.rawtrial    = 'yes';
cfg.grid.filter = tmp.avg.filter;

source = ft_sourceanalysis(cfg, data);

end

function source = sourceanalysisDICS_refdip(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.dics.refdip       = [2 5 9];
cfg.keepcsd      = 'yes';
cfg.feedback     = 'none';

source = ft_sourceanalysis(cfg, data);

end

function source = sourceanalysisDICS_refchan(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.refchan      = 'MRF11';
cfg.keepcsd      = 'yes';
cfg.feedback     = 'none';

source = ft_sourceanalysis(cfg, data);

end

function source = sourceanalysisDICS_realfilter(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.realfilter   = 'yes';
cfg.keepfilter   = 'yes';
cfg.feedback     = 'none';

source = ft_sourceanalysis(cfg, data);

end

function source = sourceanalysisDICS_fixedori(data, grid, vol)

% --- HISTORICAL --- attempt forward compatibility with function handles
if ~exist('ft_sourceanalysis', 'file') && exist('sourceanalysis', 'file')
  eval('ft_sourceanalysis = @sourceanalysis;');
end

% cfg options
cfg              = [];
cfg.method       = 'dics';
cfg.vol          = vol;
cfg.channel      = 'MEG';
cfg.grid         = grid;
cfg.frequency    = 10;
cfg.fixedori     = 'yes';
cfg.feedback     = 'none';

source = ft_sourceanalysis(cfg, data);

end

