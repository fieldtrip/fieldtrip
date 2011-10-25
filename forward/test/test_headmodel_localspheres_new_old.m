function test_headmodel_localspheres_new_old

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the data which is needed

% read in the gradiometer information
cd('/home/common/matlab/fieldtrip/data/');
hdr  = ft_read_header('Subject01.ds');
grad = hdr.grad;

% specify the file for the headshape
hdmfile  = '/home/common/matlab/fieldtrip/data/Subject01.shape';

% read in the headshape
shape = ft_read_headshape(hdmfile);

%%%%%%%%%%%%%%%%%%%%%
% do the computations

% new way - high level 
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.singlesphere = 'yes';
cfg.hdmfile = hdmfile;
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
cfg.geom    = shape;
vol{4}      = ft_prepare_headmodel(cfg);

% compare the volume conductor structures
comb = nchoosek(1:numel(vol),2);

for i=1:size(comb,1)
  chk = comb(i,:);
  try
    if ~isequal(vol{chk(1)},vol{chk(2)})
      str = sprintf('combination %d %d not successful\n',chk(1),chk(2));
      error(str)
    end
  catch me
    fprintf(me.message)
  end
end

clear vol

% new way - high level 
cfg         = [];
cfg.method  = 'localspheres';
cfg.grad    = grad;
cfg.hdmfile = hdmfile;
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
cfg.geom    = shape;
vol{4}      = ft_prepare_headmodel(cfg);
vol{4} = rmfield(vol{4},'unit');

% compare the volume conductor structures
comb = nchoosek(1:numel(vol),2);

for i=1:size(comb,1)
  chk = comb(i,:);
  try
    if ~isequal(vol{chk(1)},vol{chk(2)})
      str = sprintf('combination %d %d not successful\n',chk(1),chk(2));
      error(str)
    end
  catch me
    fprintf(me.message)
  end
end

