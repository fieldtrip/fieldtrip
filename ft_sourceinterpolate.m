function [interp] = ft_sourceinterpolate(cfg, functional, anatomical);

% FT_SOURCEINTERPOLATE reslices and interpolates a source reconstruction
% or a statistical distribution as an overlay onto an anatomical MRI.
%
% The source volume and the anatomical volume should be expressed in the
% same coordinate sytem, i.e. either both in CTF coordinates (NAS/LPA/RPA)
% or both in SPM coordinates (AC/PC). The output volume will contain a
% resliced source and anatomical volume that can be plotted together with
% FT_SOURCEPLOT or FT_SLICEINTERP, or that can be written to file using FT_SOURCEWRITE.
%
% Use as
%   [interp] = ft_sourceinterpolate(cfg, source, mri)   or
%   [interp] = ft_sourceinterpolate(cfg, stat, mri)
% where
%   source is the output of FT_SOURCEANALYSIS
%   stat   is the output of FT_SOURCESTATISTICS
%   mri    is the output of FT_READ_FCDC_MRI or the filename of a MRI
% and cfg is a structure with any of the following fields
%   cfg.parameter     = string, default is 'all'
%   cfg.interpmethod  = 'linear', 'cubic', 'nearest' or 'spline'
%   cfg.downsample    = integer number (default = 1, i.e. no downsampling)
%
% See also FT_SOURCEANALYSIS, FT_SOURCESTATISTICS, FT_READ_FCDC_MRI

% Undocumented options
%   cfg.voxelcoord = 'yes' (default) or 'no' determines whether the
%   downsampled output anatomical MRI will have the x/y/zgrid converted or
%   the homogenous transformation matrix

% Copyright (C) 2003-2007, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

%% checkdata see below!!! %%

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'trackconfig', 'on');
cfg = checkconfig(cfg, 'unused',  {'keepinside'});

% set the defaults
if ~isfield(cfg, 'parameter'),    cfg.parameter    = 'all';     end
if ~isfield(cfg, 'interpmethod'); cfg.interpmethod = 'linear';  end
if ~isfield(cfg, 'downsample');   cfg.downsample   = 1;         end
if ~isfield(cfg, 'voxelcoord'),   cfg.voxelcoord   = 'yes';     end
if ~isfield(cfg, 'feedback'),     cfg.feedback     = 'text';    end
% if ~isfield(cfg, 'sourceunits');  cfg.sourceunits  = [];        end % this is deprecated, since now autodetermined
% if ~isfield(cfg, 'mriunits');     cfg.mriunits     = [];        end % this is deprecated, since now autodetermined
cfg = checkconfig(cfg, 'deprecated', {'sourceunits', 'mriunits'});

if ischar(anatomical)
  % read the anatomical MRI data from file
  fprintf('reading MRI from file\n');
  anatomical = read_mri(anatomical);
end

% check if the input data is valid for this function and ensure that the structures correctly describes a volume
functional = checkdata(functional, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunits', 'yes');
anatomical = checkdata(anatomical, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunits', 'yes');

if isfield(cfg, 'sourceunits') && ~isempty(cfg.sourceunits)
  % this uses a deprecated option
  if ~strcmp(functional.unit, cfg.sourceunits)
    warning('the automatically determined sourceunits (%s) do not match your specification (%s)', functional.unit, cfg.sourceunits);
    functional.unit = cfg.sourceunits; % override the automatically determined units
  end
end

if isfield(cfg, 'mriunits') && ~isempty(cfg.mriunits)
  % this uses a deprecated option
  if ~strcmp(anatomical.unit, cfg.mriunits)
    warning('the automatically determined mriunits (%s) do not match your specification (%s)', anatomical.unit, cfg.sourceunits);
    anatomical.unit = cfg.mriunits; % override the automatically determined units
  end
end

if ~strcmp(functional.unit, anatomical.unit)
  fprintf('converting functional data from %s into %s\n', functional.unit, anatomical.unit);
  functional = convert_units(functional, anatomical.unit);
end

% select the parameters that should be interpolated
cfg.parameter = parameterselection(cfg.parameter, functional);
cfg.parameter = setdiff(cfg.parameter, 'inside'); % inside is handled seperately

if ~isequal(cfg.downsample, 1)
  % downsample the anatomical volume
  tmpcfg = [];
  tmpcfg.downsample = cfg.downsample;
  tmpcfg.parameter  = 'anatomy';
  anatomical = volumedownsample(tmpcfg, anatomical);
end

% collect the functional volumes that should be converted
vol_name = {};
vol_data = {};
for i=1:length(cfg.parameter)
  if ~iscell(getsubfield(functional, cfg.parameter{i}))
    vol_name{end+1} = cfg.parameter{i};
    vol_data{end+1} = getsubfield(functional, cfg.parameter{i});
  else
    fprintf('not interpolating %s, since it is not a scalar field\n', cfg.parameter{i});
  end
end

% remember the coordinate trandsformation for both
transform_ana = anatomical.transform;
transform_fun = functional.transform;

% convert the anatomical voxel positions into voxel indices into the functional volume
anatomical.transform = functional.transform \ anatomical.transform;
functional.transform = eye(4);

[fx, fy, fz] = voxelcoords(functional);
[ax, ay, az] = voxelcoords(anatomical);

% estimate the subvolume of the anatomy that is spanned by the functional volume
minfx = 1;
minfy = 1;
minfz = 1;
maxfx = functional.dim(1);
maxfy = functional.dim(2);
maxfz = functional.dim(3);
sel = ax(:)>=minfx & ...
  ax(:)<=maxfx & ...
  ay(:)>=minfy & ...
  ay(:)<=maxfy & ...
  az(:)>=minfz & ...
  az(:)<=maxfz;
fprintf('selecting subvolume of %.1f%%\n', 100*sum(sel)./prod(anatomical.dim));

% start with an empty output structure
interp = [];

dimf  = [functional.dim 1 1];
allav = zeros([anatomical.dim dimf(4:end)]);
functional.inside = functional.inside(:,:,:,1,1);

if all(functional.inside(:))
  % keep all voxels marked as inside
  interp.inside = true(anatomical.dim);
else
  % reslice and interpolate inside
  interp.inside = zeros(anatomical.dim);
  % interpolate with method nearest
  interp.inside( sel) = my_interpn(double(functional.inside), ax(sel), ay(sel), az(sel), 'nearest', cfg.feedback);
  interp.inside(~sel) = 0;
  interp.inside = logical(interp.inside);
end

% prepare the grid that is used in the interpolation
fg = [fx(:) fy(:) fz(:)];
clear fx fy fz

% reslice and interpolate all functional volumes
for i=1:length(vol_name)
  fprintf('reslicing and interpolating %s\n', vol_name{i});
  for k=1:dimf(4)
    for m=1:dimf(5)
      fv = vol_data{i}(:,:,:,k,m);
      if ~isa(fv, 'double')
        % only convert if needed, this saves memory
        fv = double(fv);
      end
      av = zeros(anatomical.dim);
      % av( sel) = my_interpn(fx, fy, fz, fv, ax(sel), ay(sel), az(sel), cfg.interpmethod, cfg.feedback);
      if islogical(vol_data{i})
        % interpolate always with method nearest
        av( sel) = my_interpn(fv, ax(sel), ay(sel), az(sel), 'nearest', cfg.feedback);
        av = logical(av);
      else
        if ~all(functional.inside(:))
          % extrapolate the outside of the functional volumes for better interpolation at the edges
          fv(~functional.inside) = griddatan(fg(functional.inside(:), :), fv(functional.inside(:)), fg(~functional.inside(:), :), 'nearest');
        end
        % interpolate functional onto anatomical grid
        av( sel) = my_interpn(fv, ax(sel), ay(sel), az(sel), cfg.interpmethod, cfg.feedback);
        clear fv
        av(~sel) = nan;
        av(~interp.inside) = nan;
      end
      allav(:,:,:,k,m) = av;
      clear av
    end
  end
  interp = setsubfield(interp, vol_name{i}, allav);
end

% add the other parameters to the output
interp.dim       = anatomical.dim;
interp.transform = transform_ana; % the original coordinate system
if ~any(strcmp(cfg.parameter, 'anatomy'))
  % copy the anatomy into the functional data
  interp.anatomy   = anatomical.anatomy;
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
cfg.previous = [];
try, cfg.previous{1} = functional.cfg; end
try, cfg.previous{2} = anatomical.cfg; end
% remember the exact configuration details in the output
interp.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this function computes the location of all voxels in head
% coordinates in a memory efficient manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, z] = voxelcoords(volume)
dim       = volume.dim;
transform = volume.transform;
if isfield(volume, 'xgrid')
  xgrid = volume.xgrid;
  ygrid = volume.ygrid;
  zgrid = volume.zgrid;
else
  xgrid = 1:dim(1);
  ygrid = 1:dim(2);
  zgrid = 1:dim(3);
end
npix = prod(dim(1:2));  % number of voxels in a single slice
nvox = prod(dim);
x = zeros(dim);
y = zeros(dim);
z = zeros(dim);
X = zeros(1,npix);
Y = zeros(1,npix);
Z = zeros(1,npix);
E = ones(1,npix);
% determine the voxel locations per slice
for i=1:dim(3)
  [X(:), Y(:), Z(:)] = ndgrid(xgrid, ygrid, zgrid(i));
  tmp = transform*[X; Y; Z; E];
  x((1:npix)+(i-1)*npix) = tmp(1,:);
  y((1:npix)+(i-1)*npix) = tmp(2,:);
  z((1:npix)+(i-1)*npix) = tmp(3,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for memory efficient interpolation
% the only reason for this function is that it does the interpolation in smaller chuncks
% this prevents memory problems that I often encountered here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [av] = my_interpn(fx, fy, fz, fv, ax, ay, az, interpmethod, feedback);
function [av] = my_interpn(fv, ax, ay, az, interpmethod, feedback);
num = numel(ax);            % total number of voxels
blocksize = floor(num/20);  % number of voxels to interpolate at once, split it into 20 chuncks
lastblock = 0;              % boolean flag for while loop
sel = 1:blocksize;          % selection of voxels that are interpolated, this is the first chunck
progress('init', feedback, 'interpolating');
while (1)
  progress(sel(1)/num, 'interpolating %.1f%%\n', 100*sel(1)/num);
  if sel(end)>num
    sel = sel(1):num;
    lastblock = 1;
  end
  av(sel) = interpn(fv, ax(sel), ay(sel), az(sel), interpmethod);
  if lastblock
    break
  end
  sel = sel + blocksize;
end
progress('close');

