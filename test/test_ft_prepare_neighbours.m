function test_ft_prepare_neighbours

% MEM 1gb
% WALLTIME 00:10:00

% TEST ft_prepare_neighbours

datainfo = ref_datasets;

% get an MEG and an EEG set (hard-coded
eeginfo = datainfo(4);
meginfo = datainfo(7);

megdir = dccnpath(meginfo.origdir);

% do the MEG processing
fname = fullfile(megdir,'latest', 'raw',meginfo.type,['preproc_' meginfo.datatype]);
data = [];
load(fname);

%% triangulation method
cfg = [];
cfg.method = 'triangulation';
neighbours = ft_prepare_neighbours(cfg, data);

% do only on some channels
cfg.channel = data.label(1:end/2);
neighbours = ft_prepare_neighbours(cfg, data);

% remove some channels
data1 = data;
data1.label{1} = [data.label{1} '_disabled'];
neighbours = ft_prepare_neighbours(cfg, data1);

%% template method
cfg = [];
cfg.method = 'template';
tic
neighbours = ft_prepare_neighbours(cfg, data);
toc

cfg = [];
cfg.method = 'template';
cfg.template = 'bti248_neighb.mat';
tic
neighbours = ft_prepare_neighbours(cfg);
toc

% do only on some channels
cfg.channel = data.label(1:end/2);
neighbours = ft_prepare_neighbours(cfg, data);

% remove some channels
data1 = data;
data1.label{1} = [data.label{1} '_disabled'];
neighbours = ft_prepare_neighbours(cfg, data1);


%% do the EEG processing
eegdir = dccnpath(eeginfo.origdir);

fname = fullfile(eegdir,'latest', 'raw',eeginfo.type,['preproc_' eeginfo.datatype]);
data = [];
load(fname);

% define elec directory
ft_pos = strfind(eeginfo.origdir, 'fieldtrip') + numel('fieldtrip');
elec_dir = fullfile(eeginfo.origdir(1:ft_pos), 'template', 'electrode');

%% triangulation method
cfg = [];
cfg.elecfile = fullfile(elec_dir, 'standard_1005.elc');
cfg.method = 'triangulation';
cfg.feedback = 'yes';
neighbours = ft_prepare_neighbours(cfg, data);

% do only on some channels
cfg.channel = data.label(1:end/2);
neighbours = ft_prepare_neighbours(cfg, data);

% remove some channels
data1 = data;
data1.label{1} = [data.label{1} '_disabled'];
neighbours = ft_prepare_neighbours(cfg, data1);

%% template method
try
  cfg = [];
  cfg.method = 'template';
  tic
  neighbours = ft_prepare_neighbours(cfg, data);
  toc
  error('no template defined but no error was thrown...')
end


cfg = [];
cfg.method = 'template';
cfg.template = 'elec1005_neighb.mat';
tic
neighbours = ft_prepare_neighbours(cfg);
toc

% do only on some channels
cfg.channel = data.label(1:end/2);
neighbours = ft_prepare_neighbours(cfg);

try
  % with data there should be an error, cause template cannot define
  % neighbours for a subset of its channels
  neighbours = ft_prepare_neighbours(cfg, data);
  error('No error thrown although there should have been one (there was one in the past if no neighbours could be defined)');
catch e
  if strfind(lower(e.message), 'no neighbours') < 0
    error('Error that was thrown not appropriately named, should have contained something with ''no neighbours'' (found/defined/whatever)');
  end
end

% use an appropriate template and remove some channels
data1 = data;
data1.label{1} = [data.label{1} '_disabled'];
neighbours = ft_prepare_neighbours(cfg, data1);

%  use an appropriate template and do only on some channels
cfg = [];
cfg.method = 'template';
cfg.template = 'elec1020_neighb.mat';
cfg.channel = data.label(1:end/2);
neighbours = ft_prepare_neighbours(cfg, data);
neighbours = ft_prepare_neighbours(cfg, data1);

