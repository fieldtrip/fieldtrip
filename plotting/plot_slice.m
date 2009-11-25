function plot_slice(data,varargin)
%
% PLOT_SLICE visualizes the slices of a MRI 3D volume
%
% Use as
%   plot_slice(data, varargin)
%
% Some defaults for the additional arguments:
%
%   'nslices'       = number of slices, (default = 4)
%   'slicerange'    = range of slices in data, (default = 'auto')
%                       'auto', full range of data
%                       [min max], coordinates of first and last slice in voxels
%   'slicedim'      = dimension to slice 1 (x-axis) 2(y-axis) 3(z-axis) (default = 3)
%   'sliceindex'    = progressive integer index of the slice (1..N) in the selected dimension
%   'title'         = string, title of the figure window
%   'colorbar'      = 'yes' or 'no' (default = 'yes')
%   'map'           = colormap assigned to the slices (default='gray')
%   'transform'     = transformation matrix from voxels to mm (default = eye(4))
%   'flat2D'        = flat multi-slice representation (default=false)
%   'tag'
%   'funcolorlim'
%   'funcolormap'
%   'opacitylim'
%   'opacitymap'
%
% Example
%   mri = read_mri('Subject01.mri');
%   figure, plot_slice(mri,'title','3D volume','colorbar','yes','map','jet')
%
% Copyright (C) 2009, Cristiano Micheli 
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get the optional input arguments
slicerange    = keyval('slicerange',  varargin); if isempty(slicerange),slicerange='auto'; end
nslices       = keyval('nslices',  varargin);  if isempty(nslices),nslices=4; end
slicedim      = keyval('slicedim',  varargin);  if isempty(slicedim),slicedim=3; end
sliceindex    = keyval('sliceindex',  varargin); if isempty(sliceindex),sliceindex=[]; end
title_        = keyval('title',  varargin);  if isempty(title_),title_=''; end
colorbar1     = keyval('colorbar',  varargin);  if isempty(colorbar1),colorbar1='no'; end
funparameter  = keyval('funparameter',  varargin);  if isempty(funparameter),funparameter=[]; end
anaparameter  = keyval('anaparameter',  varargin);  if isempty(anaparameter),anaparameter='anatomy'; end
maskparameter = keyval('maskparameter',  varargin);  if isempty(maskparameter),maskparameter=[]; end
flat2D        = keyval('flat2D',  varargin); if isempty(flat2D),flat2D=false; end
map           = keyval('map',  varargin); if isempty(map),map='gray'; end
transform     = keyval('transform',  varargin); if isempty(transform),transform=eye(4); end
tag           = keyval('tag',   varargin); if isempty(tag),tag=[]; end
funcolorlim   = keyval('funcolorlim',   varargin); if isempty(funcolorlim),funcolorlim=[]; end
funcolormap   = keyval('funcolormap',   varargin); if isempty(funcolormap),funcolormap=[]; end
opacitylim    = keyval('opacitylim',   varargin); if isempty(opacitylim),opacitylim=[]; end
opacitymap    = keyval('opacitymap',   varargin); if isempty(opacitymap),opacitymap=[]; end

%%% funparameter
% has fun?
if ~isempty(funparameter)
  if issubfield(data, funparameter)
    hasfun = 1;
    fun = getsubfield(data, funparameter);
  else
    error('funparameter not found in data');
  end
else
  hasfun = 0;
  fprintf('no functional parameter\n');
  fun = [];
end
if hasfun
  handle_fun(fun,funcolorlim,funcolormap)
end
%%% anaparameter
if isequal(anaparameter,'anatomy')
  if isfield(data, 'anatomy')
    hasana = 1;
    mri8  = isa(data.anatomy, 'uint8');
    mri16 = isa(data.anatomy, 'uint16');
    % convert integers to single precision float if neccessary
    if mri8 || mri16
      fprintf('converting anatomy to double\n');
      ana = double(data.anatomy);
    else
      ana = data.anatomy;
    end
  else
    warning('no anatomical volume present, not plotting anatomy\n')
    hasana = 0;
  end
elseif isempty(anaparameter);
  hasana = 0;
  fprintf('not plotting anatomy\n');
else
  warning('do not understand anaparameter, not plotting anatomy\n')
  hasana = 0;
end
%%% maskparameter
% has mask?
if ~isempty(maskparameter)
  if issubfield(data, maskparameter)
    if ~hasfun
      error('you can not have a mask without functional data')
    else
      hasmsk = 1;
      msk = getsubfield(data, maskparameter);
      if islogical(msk) %otherwise sign() not posible
        msk = double(msk);
      end
    end
  else
    error('maskparameter not found in data');
  end
else
  hasmsk = 0;
  fprintf('no masking parameter\n');
  msk = [];
end
if hasmsk
  handle_msk(msk,opacitylim,funparameter,maskparameter,funcolorlim,opacitymap)
end

 
%%%%% select slices
ss = setdiff([1 2 3], slicedim);

if ~isstr(slicerange)
  ind_fslice = slicerange(1);
  ind_lslice = slicerange(2);
elseif isequal(slicerange, 'auto')
  if hasfun %default
    if isfield(data,'inside')
      ind_fslice = min(find(max(max(data.inside,[],ss(1)),[],ss(2))));
      ind_lslice = max(find(max(max(data.inside,[],ss(1)),[],ss(2))));
    else
      ind_fslice = min(find(~isnan(max(max(fun,[],ss(1)),[],ss(2)))));
      ind_lslice = max(find(~isnan(max(max(fun,[],ss(1)),[],ss(2)))));
    end
  elseif hasana %if only ana, no fun
    ind_fslice = min(find(max(max(ana,[],ss(1)),[],ss(2))));
    ind_lslice = max(find(max(max(ana,[],ss(1)),[],ss(2))));
  else
    error('no functional parameter and no anatomical parameter, can not plot');
  end
else
  error('do not understand slicerange');
end

if nslices==1
  ind_allslice = (ind_fslice+ind_lslice)/2;
elseif nslices==2
  ind_allslice = [ind_fslice+(ind_lslice-ind_fslice)/5 ind_fslice+4*(ind_lslice-ind_fslice)/5];
else
  ind_allslice = linspace(ind_fslice,ind_lslice,nslices);
end
ind_allslice = round(ind_allslice);


% if i want to plot a 2D representation of several slices
if flat2D 
    plot_flatslice(hasana,hasfun,hasmsk,ana,fun,msk,slicedim,ind_allslice,colorbar1)
else
  if isempty(sliceindex)
    plot_slice_sub(ana,fun,msk,slicedim,ind_allslice,map,transform); 
  else
    plot_slice_sub(ana,fun,msk,slicedim,sliceindex,map,transform);
  end 
end
if ~isempty(title_), title(title_); end



function plot_slice_sub(data,fun,msk,slicedim,ind_allslice,map,transform) 
if ~ishold, hold on, end
ds = size(data); 

% determine location of each anatomical voxel in its own voxel coordinates
i = 1:ds(1);
j = 1:ds(2);
k = 1:ds(3);
[I, J, K] = ndgrid(i, j, k);
ijk = [I(:) J(:) K(:) ones(prod(ds),1)]';

% determine location of each anatomical voxel in head coordinates
xyz = transform * ijk;
% xyz = permute(xyz,[2 1 3]);
% xdata = reshape(xyz(:,1), [ds(2) ds(1) ds(3)]);
% ydata = reshape(xyz(:,2), [ds(2) ds(1) ds(3)]);
% zdata = reshape(xyz(:,3), [ds(2) ds(1) ds(3)]);
xdata = reshape(xyz(1,:), [ds(1) ds(2) ds(3)]);
ydata = reshape(xyz(2,:), [ds(1) ds(2) ds(3)]);
zdata = reshape(xyz(3,:), [ds(1) ds(2) ds(3)]);


if slicedim == 3 
  for i=1:length(ind_allslice)
    cdata  = squeeze(data(:,:,ind_allslice(i)));
    xdata_ = squeeze(xdata(:,:,ind_allslice(i)));
    ydata_ = squeeze(ydata(:,:,ind_allslice(i)));
    zdata_ = squeeze(zdata(:,:,ind_allslice(i)));
    news   = surface('cdata',cdata,'alphadata',cdata, 'xdata',xdata_, 'ydata',ydata_, 'zdata',zdata_);  
    set(news,'facec','interp','edgec','n','facea',0.5); 
  end
elseif slicedim == 2
  for i=1:length(ind_allslice)
    cdata  = squeeze(data(:,ind_allslice(i),:));
    xdata_ = squeeze(xdata(:,ind_allslice(i),:));
    ydata_ = squeeze(ydata(:,ind_allslice(i),:));
    zdata_ = squeeze(zdata(:,ind_allslice(i),:));
    news   = surface('cdata',cdata,'alphadata',cdata, 'xdata',xdata_, 'ydata',ydata_, 'zdata',zdata_);  
    set(news,'facec','interp','edgec','n','facea',0.5); 
  end    
elseif slicedim == 1
  for i=1:length(ind_allslice)
    cdata  = squeeze(data(ind_allslice(i),:,:));
    xdata_ = squeeze(xdata(ind_allslice(i),:,:));
    ydata_ = squeeze(ydata(ind_allslice(i),:,:));
    zdata_ = squeeze(zdata(ind_allslice(i),:,:));
    news   = surface('cdata',cdata,'alphadata',cdata, 'xdata',xdata_, 'ydata',ydata_, 'zdata',zdata_);  
    set(news,'facec','interp','edgec','n','facea',0.5); 
  end    
end
view(45,45)
colormap(map)
axis off
axis vis3d
axis equal

function plot2D(vols2D, scales)
cla;
% put 2D volumes in fun, ana and msk
hasana = length(vols2D)>0 && ~isempty(vols2D{1});
hasfun = length(vols2D)>1 && ~isempty(vols2D{2});
hasmsk = length(vols2D)>2 && ~isempty(vols2D{3});

% the transpose is needed for displaying the matrix using the Matlab image() function
if hasana; ana = vols2D{1}'; end;
if hasfun; fun = vols2D{2}'; end;
if hasmsk; msk = vols2D{3}'; end;


if hasana
  % scale anatomy between 0 and 1
  fprintf('scaling anatomy\n');
  amin = min(ana(:));
  amax = max(ana(:));
  ana = (ana-amin)./(amax-amin);
  clear amin amax;
  % convert anatomy into RGB values
  ana = cat(3, ana, ana, ana);
%   ha = imagesc(ana);
  plot_matrix(ana)
end
hold on

if hasfun
%   hf = imagesc(fun);
  plot_matrix(fun)
  caxis(scales{2});
  % apply the opacity mask to the functional data
  if hasmsk
    % set the opacity
    set(hf, 'AlphaData', msk)
    set(hf, 'AlphaDataMapping', 'scaled')
    alim(scales{3});
  elseif hasana
    set(hf, 'AlphaData', 0.5)
  end
end

axis equal
axis tight
axis xy

function handle_fun(fun,funcolorlim,funcolormap)
  % determine scaling min and max (fcolmin fcolmax) and funcolormap
  funmin = min(fun(:));
  funmax = max(fun(:));
  % smart lims: make from auto other string
  if isequal(funcolorlim,'auto') 
    if sign(funmin)>-1 && sign(funmax)>-1
      funcolorlim = 'zeromax';
    elseif sign(funmin)<1 && sign(funmax)<1
      funcolorlim = 'minzero';
    else
      funcolorlim = 'maxabs';
    end
  end
  if ischar(funcolorlim)
    % limits are given as string
    if isequal(funcolorlim,'maxabs')
      fcolmin = -max(abs([funmin,funmax]));
      fcolmax =  max(abs([funmin,funmax]));
      if isequal(funcolormap,'auto'); funcolormap = 'jet'; end;
    elseif isequal(funcolorlim,'zeromax')
      fcolmin = 0;
      fcolmax = funmax;
      if isequal(funcolormap,'auto'); funcolormap = 'hot'; end;
    elseif isequal(funcolorlim,'minzero')
      fcolmin = funmin;
      fcolmax = 0;
      if isequal(funcolormap,'auto'); funcolormap = 'cool'; end;
    else
      error('do not understand cfg.funcolorlim');
    end
  else
    if ~isempty(funcolorlim)
      % limits are numeric
      fcolmin = funcolorlim(1);
      fcolmax = funcolorlim(2);
    end
    % smart colormap
    if isequal(funcolormap,'auto') 
      if sign(fcolmin) == -1 && sign(fcolmax) == 1
        funcolormap = 'jet';
      else
        if fcolmin < 0
          funcolormap = 'cool';
        else
          funcolormap = 'hot';
        end
      end
    end
  end %if ischar
  clear funmin funmax;
  % ensure that the functional data is real
  if ~isreal(fun)
    fprintf('taking absolute value of complex data\n');
    fun = abs(fun);
  end

function handle_msk(msk,opacitylim,funparameter,maskparameter,funcolorlim,opacitymap)
  mskmin = min(msk(:));
  mskmax = max(msk(:));
  % determine the opacity limits and the opacity map
  % smart lims: make from auto other string, or equal to funcolorlim if funparameter == maskparameter
  if isequal(opacitylim,'auto')
    if isequal(funparameter,maskparameter)
      opacitylim = funcolorlim;
    else
      if sign(mskmin)>-1 && sign(mskmax)>-1
        opacitylim = 'zeromax';
      elseif sign(mskmin)<1 && sign(mskmax)<1
        opacitylim = 'minzero';
      else
        opacitylim = 'maxabs';
      end
    end
  end
  if ischar(opacitylim)
    % limits are given as string
    switch opacitylim
      case 'zeromax'
        opacmin = 0;
        opacmax = mskmax;
        if isequal(opacitymap,'auto'), opacitymap = 'rampup'; end;
      case 'minzero'
        opacmin = mskmin;
        opacmax = 0;
        if isequal(opacitymap,'auto'), opacitymap = 'rampdown'; end;
      case 'maxabs'
        opacmin = -max(abs([mskmin, mskmax]));
        opacmax =  max(abs([mskmin, mskmax]));
        if isequal(opacitymap,'auto'), opacitymap = 'vdown'; end;
      otherwise
        error('incorrect specification of cfg.opacitylim');
    end
  else
    if ~isempty(opacitylim)
      % limits are numeric
      opacmin = opacitylim(1);
      opacmax = opacitylim(2);
    end
    if isequal(opacitymap,'auto')
      if sign(opacmin)>-1 && sign(opacmax)>-1
        opacitymap = 'rampup';
      elseif sign(opacmin)<1 && sign(opacmax)<1
        opacitymap = 'rampdown';
      else
        opacitymap = 'vdown';
      end
    end
  end % handling opacitylim and opacitymap
  clear mskmin mskmax;
    
function plot_flatslice(hasana,hasfun,hasmsk,ana,fun,msk,slicedim,ind_allslice,colorbar1)
    
% make new ana, fun, msk, mskana with only the slices that will be plotted (slice dim is always third dimension)
if slicedim == 3
  if hasana; new_ana = ana(:,:,ind_allslice); clear ana; ana=new_ana; clear new_ana; end;
  if hasfun; new_fun = fun(:,:,ind_allslice); clear fun; fun=new_fun; clear new_fun; end;
  if hasmsk; new_msk = msk(:,:,ind_allslice); clear msk; msk=new_msk; clear new_msk; end;
elseif slicedim == 2
  if hasana; new_ana = ana(:,ind_allslice,:); clear ana; ana=new_ana; clear new_ana; end;
  if hasfun; new_fun = fun(:,ind_allslice,:); clear fun; fun=new_fun; clear new_fun; end;
  if hasmsk; new_msk = msk(:,ind_allslice,:); clear msk; msk=new_msk; clear new_msk; end;
elseif slicedim == 1
  if hasana; new_ana = ana(ind_allslice,:,:); clear ana; ana=new_ana; clear new_ana; end;
  if hasfun; new_fun = fun(ind_allslice,:,:); clear fun; fun=new_fun; clear new_fun; end;
  if hasmsk; new_msk = msk(ind_allslice,:,:); clear msk; msk=new_msk; clear new_msk; end;
else
  error('Error: incorrect slice dimension specification')
end
%if hasmskana; new_mskana = mskana(:,:,ind_allslice); clear mskana; mskana=new_mskana; clear new_mskana; end;

% update the dimensions of the volume
if hasana; dim=size(ana); else dim=size(fun); end;

%%%%% make "quilts", that contain all slices on 2D patched sheet
% Number of patches along sides of Quilt (M and N)
% Size (in voxels) of side of patches of Quilt (m and n)

if slicedim == 3
  m = dim(1);
  n = dim(2);
  if length(dim)==2
    dim(3) = 1;
  end
  M = ceil(sqrt(dim(3)));
  N = ceil(sqrt(dim(3)));
elseif slicedim == 2
  m = dim(1);
  n = dim(3);
  M = ceil(sqrt(dim(2)));
  N = ceil(sqrt(dim(2)));
elseif slicedim == 1
  m = dim(2);
  n = dim(3);
  M = ceil(sqrt(dim(1)));
  N = ceil(sqrt(dim(1)));
end

num_patch = N*M;
%   if slicedim~=3
%     error('only supported for slicedim=3');
%   end
num_slice = (dim(slicedim));
num_empt = num_patch-num_slice;

% put empty slides on ana, fun, msk, mskana to fill Quilt up
if slicedim == 3
  if hasana; ana(:,:,end+1:num_patch)=0; end;
  if hasfun; fun(:,:,end+1:num_patch)=0; end;
  if hasmsk; msk(:,:,end+1:num_patch)=0; end;
elseif slicedim == 2
  if hasana; ana(:,end+1:num_patch,:)=0; end;
  if hasfun; fun(:,end+1:num_patch,:)=0; end;
  if hasmsk; msk(:,end+1:num_patch,:)=0; end;
elseif slicedim == 1
  if hasana; ana(end+1:num_patch,:,:)=0; end;
  if hasfun; fun(end+1:num_patch,:,:)=0; end;
  if hasmsk; msk(end+1:num_patch,:,:)=0; end;
end

%if hasmskana; mskana(:,:,end:num_patch)=0; end;
% put the slices in the quilt
for iSlice = 1:num_slice
  xbeg = floor((iSlice-1)./M);
  ybeg = mod(iSlice-1, M);
  if slicedim == 3
    if hasana
      quilt_ana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(ana(:,:,iSlice));
    end
    if hasfun
      quilt_fun(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(fun(:,:,iSlice));
    end
    if hasmsk
      quilt_msk(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(msk(:,:,iSlice));
    end
  elseif slicedim == 2
    if hasana
      quilt_ana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(ana(:,iSlice,:));
    end
    if hasfun
      quilt_fun(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(fun(:,iSlice,:));
    end
    if hasmsk
      quilt_msk(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(msk(:,iSlice,:));
    end
  elseif slicedim == 1
    if hasana
      quilt_ana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(ana(iSlice,:,:));
    end
    if hasfun
      quilt_fun(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(fun(iSlice,:,:));
    end
    if hasmsk
      quilt_msk(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(msk(iSlice,:,:));
    end
  end

  %     if hasmskana
  %       quilt_mskana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(mskana(:,:,iSlice));
  %     end
end
% make vols and scales, containes volumes to be plotted (fun, ana, msk) %added ingnie
if hasana; vols2D{1} = quilt_ana; scales{1} = []; end; % needed when only plotting ana
if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end;
if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end;

plot2D(vols2D, scales);
axis off

if strcmp(colorbar1,  'yes'),
  if hasfun
    % use a normal Matlab colorbar
    hc = colorbar;
    set(hc, 'YLim', [fcolmin fcolmax]);
  else
    warning('no colorbar possible without functional data')
  end
end
    