function test_issue1067

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_prepare_layout ft_prepare_mesh ft_prepare_sourcemodel ft_prepare_headmodel ft_prepare_neighbours

% For some fields there are multiple ways to obtain them. The desired use in the
% high-level functions is to call the function that is the highest in the hierarchy.
% 
% ft_prepare_xxx -> ft_fetch_xxx -> ft_read_xxx/loadvar

%%

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1067'));


%% layout

cfg = [];
cfg.layout = 'CTF151.lay';
cfg.layout = ft_prepare_layout(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.layout));
cfg.layout = ft_prepare_layout(cfg); % using FT_PREPARE_XXX again

cfg = [];
cfg.layout = 'CTF151.mat';
cfg.layout = ft_prepare_layout(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.layout));
cfg.layout = ft_prepare_layout(cfg); % using FT_PREPARE_XXX again


%% headshape

cfg = [];
cfg.headshape = 'Subject01.shape';
cfg.headshape = ft_prepare_mesh(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.headshape));
cfg.headshape = ft_prepare_mesh(cfg); % using FT_PREPARE_XXX again

cfg = [];
cfg.headshape = 'Subject01_headshape.mat';
cfg.headshape = ft_prepare_mesh(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.headshape));
cfg.headshape = ft_prepare_mesh(cfg); % using FT_PREPARE_XXX again


%% sourcemodel

cfg = [];
cfg.sourcemodel = 'Subject01_sourcemodel.mat';
cfg.sourcemodel = ft_prepare_sourcemodel(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.sourcemodel));
cfg.sourcemodel = ft_prepare_sourcemodel(cfg); % using FT_PREPARE_XXX again


%% headmodel

cfg = [];
cfg.headmodel = 'Subject01.hdm';
cfg.headmodel = ft_prepare_headmodel(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.headmodel));
cfg.headmodel = ft_prepare_headmodel(cfg); % using FT_PREPARE_XXX again

cfg = [];
cfg.headmodel = 'Subject01_headmodel.mat';
cfg.headmodel = ft_prepare_headmodel(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.headmodel));
cfg.headmodel = ft_prepare_headmodel(cfg); % using FT_PREPARE_XXX again


%% neighbours

cfg = [];
cfg.neighbours = 'ctf151_neighb.mat';
cfg.neighbours = ft_prepare_neighbours(cfg); % using FT_PREPARE_XXX
assert(isstruct(cfg.neighbours));
cfg.neighbours = ft_prepare_neighbours(cfg); % using FT_PREPARE_XXX again


%% grad

cfg = [];
cfg.grad = 'Subject01_grad.mat';
cfg.grad = ft_fetch_sens(cfg); % using FT_FETCH_XXX
assert(isstruct(cfg.grad));
cfg.grad = ft_fetch_sens(cfg); % using FT_FETCH_XXX again


%% elec

cfg = [];
cfg.elec = 'Subject01_elec.mat';
cfg.elec = ft_fetch_sens(cfg); % using FT_FETCH_XXX
assert(isstruct(cfg.elec));
cfg.elec = ft_fetch_sens(cfg); % using FT_FETCH_XXX again


%% opto

cfg = [];
cfg.opto = 'opto.mat';
cfg.opto = ft_fetch_sens(cfg); % using FT_FETCH_XXX
assert(isstruct(cfg.opto));
cfg.opto = ft_fetch_sens(cfg); % using FT_FETCH_XXX again


%% event

cfg = [];
cfg.event = 'Subject01_event.mat';
cfg.event = ft_read_event(cfg.event); % using FT_READ_XXX
assert(isstruct(cfg.event));


%% mri

cfg = [];
cfg.mri = 'Subject01.mri';
cfg.mri = ft_read_mri(cfg.mri);  % using FT_READ_XXX
assert(isstruct(cfg.mri));

cfg = [];
cfg.mri = 'Subject01_mri.mat';
cfg.mri = ft_read_mri(cfg.mri);  % using FT_READ_XXX
assert(isstruct(cfg.mri));



