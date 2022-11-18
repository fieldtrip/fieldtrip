function [resliced] = ft_volumereslice(cfg, mri)

% FT_VOLUMERESLICE flips, permutes, interpolates and reslices a volume along the
% principal axes of the coordinate system according to a specified resolution.
%
% Use as
%   mri = ft_volumereslice(cfg, mri)
% where the input MRI should be a single anatomical or functional MRI volume that
% results from FT_READ_MRI or FT_VOLUMEREALIGN. You can visualize the the input and
% output using FT_SOURCEPLOT.
%
% The configuration structure can contain
%   cfg.method     = string, 'flip', 'nearest', 'linear', 'cubic' or 'spline' (default = 'linear')
%   cfg.downsample = integer number (default = 1, i.e. no downsampling)
%
% If you specify the method as 'flip', it will only permute and flip the volume, but
% not perform any interpolation. For the other methods the input volumetric data will
% also be interpolated on a regular voxel grid.
%
% For the interpolation methods you should specify
%   cfg.resolution = number, in units of distance (e.g. mm)
%   cfg.xrange     = [min max], in units of distance (e.g. mm)
%   cfg.yrange     = [min max], in units of distance (e.g. mm)
%   cfg.zrange     = [min max], in units of distance (e.g. mm)
% or alternatively with
%   cfg.dim        = [nx ny nz], size of the volume in each direction
%
% If the input MRI has a coordsys-field and you don't specify explicit the
% xrange/yrange/zrange, the centre of the volume will be shifted (with respect to the
% origin of the coordinate system), for the brain to fit nicely in the box.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_VOLUMEREALIGN, FT_VOLUMEDOWNSAMPLE, FT_SOURCEINTERPOLATE, FT_SOURCEPLOT

% Copyright (C) 2010-2020, Robert Oostenveld & Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    mri
ft_preamble provenance mri

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function and ensure that the structures correctly describes a volume
if isfield(mri, 'inside')
  mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes', 'insidestyle', 'logical');
else
  mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes');
end

% set the defaults
cfg.method     = ft_getopt(cfg, 'method', 'linear');
cfg.downsample = ft_getopt(cfg, 'downsample', 1);

if isequal(cfg.method, 'flip')
  % these do not apply when flipping
  cfg = ft_checkconfig(cfg, 'forbidden', {'resolution', 'xrange', 'yrange', 'zrange', 'dim'});
else
  % these only applies when interpolating
  cfg.resolution = ft_getopt(cfg, 'resolution', 1 * ft_scalingfactor('mm', mri.unit)); % default is 1 mm, but the actual number depends on the units. See bug2906
  cfg.xrange     = ft_getopt(cfg, 'xrange', []);
  cfg.yrange     = ft_getopt(cfg, 'yrange', []);
  cfg.zrange     = ft_getopt(cfg, 'zrange', []);
  cfg.dim        = ft_getopt(cfg, 'dim', []);
  
  if isfield(mri, 'coordsys')
    % use some prior knowledge to optimize the location of the bounding box
    % with respect to the origin of the coordinate system
    switch mri.coordsys
      case {'ctf', '4d', 'bti', 'eeglab'}
        xshift = 30./cfg.resolution;
        yshift = 0;
        zshift = 40./cfg.resolution;
      case {'neuromag', 'itab'}
        xshift = 0;
        yshift = 30./cfg.resolution;
        zshift = 40./cfg.resolution;
      case {'acpc', 'spm', 'mni', 'tal'}
        ft_warning('FIXME, the bounding box needs a better default');
        xshift = 0;
        yshift = 0;
        zshift = 0;
      otherwise
        xshift = 0;
        yshift = 0;
        zshift = 0;
    end
  else % if no coordsys is present
    xshift = 0;
    yshift = 0;
    zshift = 0;
  end
  
  if ~isempty(cfg.dim)
    xrange = [-cfg.dim(1)/2+0.5 cfg.dim(1)/2-0.5] * cfg.resolution + xshift;
    yrange = [-cfg.dim(2)/2+0.5 cfg.dim(2)/2-0.5] * cfg.resolution + yshift;
    zrange = [-cfg.dim(3)/2+0.5 cfg.dim(3)/2-0.5] * cfg.resolution + zshift;
  else % if no cfg.dim is specified, use defaults
    range = [-127.5 127.5] * cfg.resolution; % 255 mm^3 bounding box, assuming human brain
    xrange = range + xshift;
    yrange = range + yshift;
    zrange = range + zshift;
  end
  
  % if ranges have not been specified by the user
  if isempty(cfg.xrange)
    cfg.xrange = xrange;
  end
  if isempty(cfg.yrange)
    cfg.yrange = yrange;
  end
  if isempty(cfg.zrange)
    cfg.zrange = zrange;
  end
  
end % if method~=fip

if cfg.downsample~=1
  % optionally downsample the anatomical and/or functional volumes
  tmpcfg = keepfields(cfg, {'downsample', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  mri = ft_volumedownsample(tmpcfg, mri);
  % restore the provenance information
  [cfg, mri] = rollback_provenance(cfg, mri);
end

% determine the fields to reslice
fn = fieldnames(mri);
fn = setdiff(fn, {'pos', 'tri', 'inside', 'outside', 'time', 'freq', 'dim', 'transform', 'unit', 'coordsys', 'cfg', 'hdr'}); % remove fields that do not represent the data
dimord = cell(size(fn));
for i=1:numel(fn)
  dimord{i} = getdimord(mri, fn{i});
end
fn = fn(strcmp(dimord, 'dim1_dim2_dim3'));

if strcmp(cfg.method, 'flip')
  % this uses some private functions that change the volumes and the transform
  resliced = volumepermute(mri); % this makes the transform approximately diagonal
  flipvec = false(1,3);
  flipvec(1) = resliced.transform(1,1)<0;
  flipvec(2) = resliced.transform(2,2)<0;
  flipvec(3) = resliced.transform(3,3)<0;
  resliced = volumeflip(resliced, flipvec); % this flips along each of the dimensions
  
else
  % compute the desired grid positions
  xgrid = cfg.xrange(1):cfg.resolution:cfg.xrange(2);
  ygrid = cfg.yrange(1):cfg.resolution:cfg.yrange(2);
  zgrid = cfg.zrange(1):cfg.resolution:cfg.zrange(2);
  
  resliced           = [];
  resliced.dim       = [length(xgrid) length(ygrid) length(zgrid)];
  resliced.transform = translate([cfg.xrange(1) cfg.yrange(1) cfg.zrange(1)]) * scale([cfg.resolution cfg.resolution cfg.resolution]) * translate([-1 -1 -1]);
  resliced.anatomy   = zeros(resliced.dim, 'int8');
  resliced.unit      = mri.unit;
  
  % these take a lot of memory
  clear xgrid ygrid zgrid
  
  % these are the same in the resliced as in the input anatomical MRI
  if isfield(mri, 'coordsys')
    resliced.coordsys = mri.coordsys;
  end
  
  fprintf('reslicing from [%d %d %d] to [%d %d %d]\n', mri.dim(1), mri.dim(2), mri.dim(3), resliced.dim(1), resliced.dim(2), resliced.dim(3));
  
  % the actual work is being done by ft_sourceinterpolate
  % this interpolates the real volume on the resolution that is defined for the resliced volume
  tmpcfg                = [];
  tmpcfg.parameter      = fn;
  tmpcfg.interpmethod   = cfg.method;
  resliced              = ft_sourceinterpolate(tmpcfg, mri, resliced);
  resliced.cfg.previous = resliced.cfg.previous{1}; % the 2nd input is a dummy variable
  cfg.method = resliced.cfg.interpmethod;           % remember the method that was used
  % restore the provenance information
  [cfg, resliced] = rollback_provenance(cfg, resliced);
  
  % remove fields that were not present in the input
  % this applies specifically to the 'inside' field that may have been added by ft_sourceinterpolate
  resliced = keepfields(resliced, fieldnames(mri));
  
  % convert any non-finite values to 0 to avoid problems later on
  for i=1:numel(fn)
    resliced.(fn{i})(~isfinite((fn{i}))) = 0;
  end
  
end % if method=flip or interpolate

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   mri
ft_postamble provenance resliced
ft_postamble history    resliced
ft_postamble savevar    resliced
