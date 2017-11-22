function test_headmodel_localspheres_new_old

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_prepare_localspheres ft_prepare_headmodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the gradiometer information
cd(dccnpath('/home/common/matlab/fieldtrip/data'));
hdr  = ft_read_header('Subject01.ds');
grad = hdr.grad;

% specify the file for the headshape
hdmfile  = dccnpath('/home/common/matlab/fieldtrip/data/Subject01.shape');

% read in the headshape
shape = ft_read_headshape(hdmfile);

%%%%%%%%%%%%%%%%%%%%%
% do the computations

% new way - high level
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.singlesphere = 'yes';
cfg.headshape = hdmfile;
vol{1}      = ft_prepare_headmodel(cfg);

% old way
cfg              = [];
cfg.singlesphere = 'yes';
cfg.headshape    = hdmfile;
cfg.grad         = grad;
cfg.singlesphere = 'yes';
vol{2}           = ft_prepare_localspheres(cfg);

% old way
cfg              = [];
cfg.headshape    = shape;
cfg.singlesphere = 'yes';
cfg.grad         = grad;
vol{3}           = ft_prepare_localspheres(cfg);

% new way - high level
cfg         = [];
cfg.method  = 'localspheres';
cfg.singlesphere = 'yes';
cfg.grad    = grad;
vol{4}      = ft_prepare_headmodel(cfg, shape);

% compare the volume conductor structures
comb = nchoosek(1:numel(vol),2);

% compute the leadfields for a comparison
[vol{1}, grad] = ft_prepare_vol_sens(vol{1}, grad);
[vol{2}, grad] = ft_prepare_vol_sens(vol{2}, grad);
[vol{3}, grad] = ft_prepare_vol_sens(vol{3}, grad);
[vol{4}, grad] = ft_prepare_vol_sens(vol{3}, grad);
lf{1} = ft_compute_leadfield([0 10 60], grad, vol{1});
lf{2} = ft_compute_leadfield([0 10 60], grad, vol{2});
lf{3} = ft_compute_leadfield([0 10 60], grad, vol{3});
lf{4} = ft_compute_leadfield([0 10 60], grad, vol{4});

% compare the leadfields in all possible combinations
comb = nchoosek(1:numel(vol),2);
for j=1:size(comb,1)
  chk = comb(j,:);
  err = norm(lf{chk(1)} - lf{chk(2)}) / norm(lf{chk(1)});
  if err>0.001
    error('combination %d %d not successful\n',chk(1),chk(2));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear vol lf

% new way - high level
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.headshape = hdmfile;
vol{1}      = ft_prepare_headmodel(cfg);
vol{1} = rmfield(vol{1},'unit');

% old way
cfg              = [];
cfg.headshape    = hdmfile;
cfg.grad         = grad;
vol{2}           = ft_prepare_localspheres(cfg);

% old way
cfg              = [];
cfg.headshape    = shape;
cfg.grad         = grad;
vol{3}           = ft_prepare_localspheres(cfg);

% new way - high level
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
vol{4}      = ft_prepare_headmodel(cfg, shape);
vol{4} = rmfield(vol{4},'unit');

% compute the leadfields for a comparison
[vol{1}, grad] = ft_prepare_vol_sens(vol{1}, grad);
[vol{2}, grad] = ft_prepare_vol_sens(vol{2}, grad);
[vol{3}, grad] = ft_prepare_vol_sens(vol{3}, grad);
[vol{4}, grad] = ft_prepare_vol_sens(vol{3}, grad);
lf{1} = ft_compute_leadfield([0 10 60], grad, vol{1});
lf{2} = ft_compute_leadfield([0 10 60], grad, vol{2});
lf{3} = ft_compute_leadfield([0 10 60], grad, vol{3});
lf{4} = ft_compute_leadfield([0 10 60], grad, vol{4});

% compare the leadfields in all possible combinations
comb = nchoosek(1:numel(vol),2);
for j=1:size(comb,1)
  chk = comb(j,:);
  err = norm(lf{chk(1)} - lf{chk(2)}) / norm(lf{chk(1)});
  if err>0.001
    error('combination %d %d not successful\n',chk(1),chk(2));
  end
end


