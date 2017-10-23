function test_bug1166

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_headmodel ft_headmodel_bem_asa 

% This function tests that the inputs for the headmodel functions are
% read-in correctly

% ASA headmodel
filename = dccnpath('/home/common/matlab/fieldtrip/template/headmodel/skin/standard_skin_1222.vol');
cfg = [];
cfg.method = 'asa';
cfg.hdmfile = filename;
vol1 = ft_prepare_headmodel(cfg);
vol2 = ft_headmodel_asa(filename);
vol1 = rmfield(vol1,'cfg');
if ~isequal(vol1,vol2)
  error('test failed!')
end

% DIPOLI headmodel
bnd = vol1.bnd;
cfg = [];
cfg.method = 'dipoli';
vol1 = ft_prepare_headmodel(cfg,bnd);
vol1 = rmfield(vol1,'cfg');
vol1 = rmfield(vol1,'unit');
vol2 = ft_headmodel_dipoli(bnd);
if ~isequal(vol1,vol2)
  error('test failed!')
end

% Next sections obsolete, as hdmfile and headshape no longer possible
% options for input use with dipoli
%
% cfg = [];
% cfg.method = 'dipoli';
% cfg.hdmfile = filename;
% vol3 = ft_prepare_headmodel(cfg);
% vol3 = rmfield(vol3,'cfg');
% vol3 = rmfield(vol3,'unit');
% if ~isequal(vol1,vol3)
%   error('test failed!')
% end
% 
% cfg = [];
% cfg.method = 'dipoli';
% cfg.headshape = filename;
% vol4 = ft_prepare_headmodel(cfg);
% vol4 = rmfield(vol4,'cfg');
% vol4 = rmfield(vol4,'unit');
% if ~isequal(vol1,vol4)
%   error('test failed!')
% end
% 
% 
