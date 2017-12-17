function test_bug2518a

% WALLTIME 00:30:00
% MEM 2500mb

% TEST ft_componentanalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for these functions I don't really know whether and how they should work
% ft_channelrepair
% ft_combineplanar
%
% for these functions it could in principle work
% ft_megplanar
% ft_scalpcurrentdensity
%
% the following functions have not been seriously considered yet
% ft_connectivityanalysis
% ft_megrealign
% ft_mvaranalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% prepare some data
cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.demean = 'yes';
cfg.channel = 'MEG';
raw = ft_preprocessing(cfg);

% this takes less than 5 minutes
cfg = [];
cfg.method = 'runica';
cfg.runica.pca = 30;
comp = ft_componentanalysis(cfg, raw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the actual testing starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.detrend = 'no';
cfg.demean = 'yes';
cfg.resamplefs = 100;
output = ft_resampledata(cfg, comp);

assert(isfield(output, 'topo'),      'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
output = ft_timelockanalysis(cfg, comp);

assert(isfield(output, 'topo'), 'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.method = 'mtmfft';
cfg.foilim = [1 30];
cfg.taper = 'hanning';
output = ft_freqanalysis(cfg, comp);

assert(isfield(output, 'topo'), 'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
output = ft_appenddata(cfg, comp, comp);

assert(isfield(output, 'topo'),      'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');

% make an alternative but different component structure
compa = comp;
compa.topo = randn(size(compa.topo));

try
  problem = false;
  output = ft_appenddata(cfg, comp, compa);
  problem = true;
end
assert(~problem, 'different component structures should not be appended');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
output = ft_channelnormalise(cfg, comp);

assert(~isfield(output, 'topo'),      'topo is present');
assert(~isfield(output, 'topolabel'), 'topolabel is present');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.offset = 1;
output = ft_redefinetrial(cfg, comp);

assert(isfield(output, 'topo'),      'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 40;
output = ft_preprocessing(cfg, comp);

assert(isfield(output, 'topo'),      'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false
  % this should not run in non-interactive mode
  cfg = [];
  output = ft_rejectvisual(cfg, comp);
  
  assert(isfield(output, 'topo'),      'topo is missing');
  assert(isfield(output, 'topolabel'), 'topolabel is missing');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg = [];
cfg.prewindow = 0.1;
cfg.postwindow = 0.1;
output = ft_interpolatenan(cfg, comp);

assert(isfield(output, 'topo'),      'topo is missing');
assert(isfield(output, 'topolabel'), 'topolabel is missing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  problem = false;
  output = ft_megrealign(cfg, comp);
  problem = true;
end
assert(~problem, 'this function should fail on component data');



