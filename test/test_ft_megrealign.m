function test_ft_megrealign

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_megrealign ft_prepare_neighbours ft_transform_geometry

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

datainfo = ref_datasets;
sel      = match_str({datainfo.datatype},{'bti148' 'bti248' 'ctf151' 'ctf275' 'itab153' 'yokogawa160'}');
datainfo = datainfo(sel);

% as of yet, the code does not explicitly test anything, it just checks
% whether all runs through smoothly
for k = 1:numel(datainfo)
  fname = fullfile(datainfo(k).origdir,'latest/raw',datainfo(k).type,['preproc_',datainfo(k).datatype]);
  load(fname);
  
  cfg = [];
  cfg.channel = 'MEG';
  cfg.demean = 'yes';
  data = ft_preprocessing(cfg, data);
  
  vol      = [];
  vol.o    = [0 0 4];
  vol.r    = 10;
  vol.unit = 'cm';
  
  % ensure units in the gradiometer array and volume conductor to be equal
  data.grad = ft_convert_units(data.grad);
  data.grad = ft_convert_units(data.grad, 'cm');
  vol       = ft_convert_units(vol, data.grad.unit);
  
  % make 2 copies
  data2 = data;
  data3 = data;
  
  data2.grad = ft_transform_geometry([[eye(3) [0 0 -2]'];[0 0 0 1]], data.grad);
  data3.grad = ft_transform_geometry([1 0 0 0;0 cos(0.1) -sin(0.1) 0;0 sin(0.1) cos(0.1) 0;0 0 0 1], data.grad);
  
  cfg = [];
  cfg.template{1} = data.grad;
  cfg.template{2} = data2.grad;
  cfg.template{3} = data3.grad;
  cfg.inwardshift = 1;
  cfg.vol         = vol;
  
  interp = ft_megrealign(cfg, data);
  interp2 = ft_megrealign(cfg, data2);
  interp3 = ft_megrealign(cfg, data3);
  
  
end

