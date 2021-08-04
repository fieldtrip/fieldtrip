function test_issue1843

% WALLTIME 00:30:00
% MEM 4gb
% DEPENDENCY ft_sourceplot ft_plot_slice


% this test function checks whether ft_sourceplot (method: ortho) produces
% a figure (without error) using cfg.maskstyle = 'colormix'

% create a source structure with anatomical data and with functional data
anat = reshape(linspace(0,1,8),[2 2 2]);
pow  = randn(2,2,2);

source.pow = pow;
source.anatomy = anat;
source.dim = [2 2 2];
source.transform = eye(4);
source.transform(1:3,4) = -2;

cfg = [];
cfg.funparameter = 'pow';
ft_sourceplot(cfg, source);

cfg.maskstyle = 'colormix';
ft_sourceplot(cfg, source);

%source.mask = source.anatomy;
source.mask = ones(2,2,2)./2;
cfg = rmfield(cfg, 'maskstyle');
cfg.maskparameter = 'mask';
cfg.opacitylim = [0 1];
ft_sourceplot(cfg, source);

cfg.maskstyle = 'colormix';
ft_sourceplot(cfg, source);

source.mask = source.anatomy;
cfg = removefields(cfg, {'maskstyle' 'opacitylim'});
cfg.maskparameter = 'mask';
ft_sourceplot(cfg, source);

cfg.maskstyle = 'colormix';
ft_sourceplot(cfg, source);
