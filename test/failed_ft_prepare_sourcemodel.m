function failed_ft_prepare_sourcemodel

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_prepare_sourcemodel
% TEST ft_prepare_sourcemodel ft_read_sens ft_read_vol 

% function to test ft_prepare_sourcemodel given configuration options (cfg), 
% a single sphere volume condution model (vol), and gradiometer information
% (sens)
%
% additionally, we try both (1) a standard CTF275 gradiometer file and 
% (2) an extended version derived from ft_headmovement with cfg.numclusters = 10 
% A. Stolk

success = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load single sphere volume conduction model
vol = ft_read_vol(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_singlesphere.mat'), 'vol');

% load gradiometer information of an exemplary subject
grad_standard = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/sens/ctf275.mat'));


% load the same gradiometer information treated with ft_headmovement (10 clusters)
grad_extended = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/sens/ctf275_10clusters.mat'));

%%%%%%%%%%%%%%%%%%%%%
% do the computations (standard)

% create config options
cfg                 = [];
cfg.symmetry        = [];
cfg.grid.resolution = 2;

[grid, cfg] = ft_prepare_sourcemodel(cfg, vol, grad_standard);

% check whether a grid could be computed
success     = success && ~isempty(grid);
if ~success
  error('ft_prepare_sourcemodel was not able make a grid');
end

% check whether there are potential dipoles inside the brain
success     = success && ~isempty(grid.inside);
if ~success
  error('ft_prepare_sourcemodel was not able to determine the inside brain');
end

% check whether the inside field is constrained to the positions inside the
% volume conductor model
success = success && numel(grid.inside)~=size(grid.pos,1);
if ~success
  error('ft_prepare_sourcemodel was not able to constrain the inside positions');
end

%%%%%%%%%%%%%%%%%%%%%
% do the computations (extended gradiometer)

% create config options
cfg                 = [];
cfg.symmetry        = [];
cfg.grid.resolution = 2;

[grid, cfg] = ft_prepare_sourcemodel(cfg, vol, grad_extended);

% check whether a grid could be computed
success     = success && ~isempty(grid);
if ~success
  error('ft_prepare_sourcemodel was not able make a grid using extended gradiometer information');
end

% check whether there are potential dipoles inside the brain
success     = success && ~isempty(grid.inside);
if ~success
  error('ft_prepare_sourcemodel was not able to determine the inside brain using extended gradiometer information');
end

%%%%%%%%%%%%%%%%%%%%%
% clean up
clear cfg
clear vol
clear grad_standard
clear grad_extended
