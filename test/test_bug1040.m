function test_bug1040

% MEM 1gb
% WALLTIME 00:03:05

% TEST: bug1040

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
load('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_singlesphere.mat', 'vol');

% load gradiometer information of an exemplary subject
load('/home/common/matlab/fieldtrip/data/test/latest/sens/ctf275.mat', 'grad');
grad_standard = grad; clear grad;

% load the same gradiometer information treated with ft_headmovement (10
% clusters)
load('/home/common/matlab/fieldtrip/data/test/latest/sens/ctf275_10clusters.mat', 'grad');
grad_extended = grad; clear grad;

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
clear cfg;
clear vol;
clear grad_standard;
clear grad_extended;
