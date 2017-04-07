function test_bug1806

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_componentanalysis ft_rejectcomponent ft_megplanar ft_combineplanar ft_megrealign ft_datatype_sens

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

ctf151_sens = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/sens/ctf151.mat'));
ctf275_sens = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/sens/ctf275.mat'));

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.ds');
cfg.trl = [1 900 0];
data = ft_preprocessing(cfg);
% the following applies since 30 October 2012
% the type will be added by ft_datatype_sens if not present
assert(strcmp(data.grad.type, 'ctf151'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% megplanar and combineplanar

cfg = [];
cfg.method = 'distance';
neighbours = ft_prepare_neighbours(cfg, data);

cfg = [];
cfg.neighbours = neighbours;
data_p  = ft_megplanar(cfg,data);
if isfield(data_p.grad, 'type')
  % it should not be ctf151 any more
  assert(strcmp(data_p.grad.type, 'ctf151_planar'));
else
  warning('gradiometer type is missing')
  assert(ft_senstype(data_p.grad, 'ctf151_planar'));
end

cfg = [];
data_pc = ft_combineplanar(cfg,data_p);
if isfield(data_pc.grad, 'type')
  % it should again be ctf151
  assert(strcmp(data_pc.grad.type, 'ctf151_planar_combined'));
else
  warning('gradiometer type is missing');
  assert(ft_senstype(data_pc.grad, 'ctf151_planar_combined'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% megrealign with another template

vol.type = 'singlesphere';
vol.unit = 'cm';
vol.r = 12;
vol.o = [0 0 4];

cfg = [];
cfg.vol = vol;
cfg.inwardshift = 0;
cfg.template = ctf275_sens;
data_r = ft_megrealign(cfg, data);
assert(strcmp(data_r.grad.type, 'ctf275'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% componentanalysis and rejectcomponent

cfg = [];
cfg.method = 'pca';
data_c = ft_componentanalysis(cfg, data);

% this is how it used to be up to 19 Feb 2013
%
% if isfield(data_c.grad, 'type')
%   % don't know what it is, but it should not be ctf151
%   assert(~strcmp(data_c.grad.type, 'ctf151'));
% else
%   warning('gradiometer type is missing');
%   assert(~ft_senstype(data_c.grad, 'ctf151'));
% end

% on 19 Feb 2013 I changed it, because of http://bugzilla.fcdonders.nl/show_bug.cgi?id=1959#c7
if isfield(data_c.grad, 'type')
  % it should still be detected as ctf151
  assert(strcmp(data_c.grad.type, 'ctf151'));
else
  % it should still be detected as ctf151
  warning('gradiometer type is missing');
  assert(ft_senstype(data_c.grad, 'ctf151'));
end


cfg = [];
cfg.component = [];
data_cb = ft_rejectcomponent(cfg, data_c);
if isfield(data_cb.grad, 'type')
  % it should again be ctf151
  assert(strcmp(data_cb.grad.type, 'ctf151'));
else
  warning('gradiometer type is missing');
  assert(ft_senstype(data_cb.grad, 'ctf151'));
end
