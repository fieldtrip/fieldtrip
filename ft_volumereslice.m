function [resliced] = ft_volumereslice(cfg, mri)

% FT_VOLUMERESLICE interpolates and reslices a volume along the
% principal axes of the coordinate system according to a specified
% resolution.
%
% Use as
%   mri = ft_volumereslice(cfg, mri)
% where the input mri should be a single anatomical or functional MRI
% volume that was for example read with FT_READ_MRI.
%
% The configuration structure can contain
%   cfg.resolution = number, in physical units
%   cfg.xrange     = [min max], in physical units
%   cfg.yrange     = [min max], in physical units
%   cfg.zrange     = [min max], in physical units
% or alternatively with
%   cfg.dim        = [nx ny nz], size of the volume in each direction
%
% If the input mri has a coordsys-field, the centre of the volume will be
% shifted (with respect to the origin of the coordinate system), for the
% brain to fit nicely in the box.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_VOLUMEDOWNSAMPLE, FT_SOURCEINTERPOLATE

% Undocumented local options:
%   cfg.downsample

% Copyright (C) 2010-2013, Robert Oostenveld & Jan-Mathijs Schoffelen
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
ft_preamble loadvar mri
ft_preamble provenance mri
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function and ensure that the structures correctly describes a volume
if isfield(mri, 'inside')
  mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes', 'inside', 'logical');
else
  mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes');
end

% set the defaults
% set voxel resolution according to the input units- see bug2906
unitcheckmm = strcmp(mri.unit,'mm');
if unitcheckmm==1;
    cfg.resolution = ft_getopt(cfg, 'resolution', 1);
end;
unitcheckcm = strcmp(mri.unit,'cm');
if unitcheckcm==1;
    cfg.resolution = ft_getopt(cfg, 'resolution', .1);
end;
unitcheckm = strcmp(mri.unit,'m');
if unitcheckm==1;
    cfg.resolution = ft_getopt(cfg, 'resolution', .001);
end;
%cfg.resolution = ft_getopt(cfg, 'resolution', 1);
cfg.downsample = ft_getopt(cfg, 'downsample', 1);
cfg.xrange     = ft_getopt(cfg, 'xrange', []);
cfg.yrange     = ft_getopt(cfg, 'yrange', []);
cfg.zrange     = ft_getopt(cfg, 'zrange', []);
cfg.dim        = ft_getopt(cfg, 'dim', []); % alternatively use ceil(mri.dim./cfg.resolution)

if isfield(mri, 'coordsys')
  % use some prior knowledge to optimize the location of the bounding box
  % with respect to the origin of the coordinate system
  switch mri.coordsys
    case {'ctf' '4d' 'bti'}
      xshift = 30./cfg.resolution;
      yshift = 0;
      zshift = 40./cfg.resolution;
    case {'itab' 'neuromag'}
      xshift = 0;
      yshift = 30./cfg.resolution;
      zshift = 40./cfg.resolution;
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

if cfg.downsample~=1
  % optionally downsample the anatomical and/or functional volumes
  tmpcfg = keepfields(cfg, {'downsample'});
  mri = ft_volumedownsample(tmpcfg, mri);
  % restore the provenance information
  [cfg, mri] = rollback_provenance(cfg, mri);
end

% compute the desired grid positions
xgrid = cfg.xrange(1):cfg.resolution:cfg.xrange(2);
ygrid = cfg.yrange(1):cfg.resolution:cfg.yrange(2);
zgrid = cfg.zrange(1):cfg.resolution:cfg.zrange(2);

resliced           = [];
resliced.dim       = [length(xgrid) length(ygrid) length(zgrid)];
resliced.transform = translate([cfg.xrange(1) cfg.yrange(1) cfg.zrange(1)]) * scale([cfg.resolution cfg.resolution cfg.resolution]) * translate([-1 -1 -1]);
resliced.anatomy   = zeros(resliced.dim, 'int8');
resliced.unit      = mri.unit;

clear xgrid ygrid zgrid

% these are the same in the resliced as in the input anatomical MRI
if isfield(mri, 'coordsys')
  resliced.coordsys = mri.coordsys;
end

fprintf('reslicing from [%d %d %d] to [%d %d %d]\n', mri.dim(1), mri.dim(2), mri.dim(3), resliced.dim(1), resliced.dim(2), resliced.dim(3));

% the actual work is being done by ft_sourceinterpolate, which interpolates the real mri volume
% on the resolution that is defined for the resliced volume
tmpcfg = [];
tmpcfg.parameter = 'anatomy';
resliced = ft_sourceinterpolate(tmpcfg, mri, resliced);

% remove fields that were not present in the input, this applies specifically to
% the 'inside' field that may have been added by ft_sourceinterpolate
resliced = rmfield(resliced, setdiff(fieldnames(resliced), fieldnames(mri)));

% convert any non-finite values to 0 to avoid problems later on
resliced.anatomy(~isfinite(resliced.anatomy)) = 0;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance resliced
ft_postamble history    resliced
ft_postamble savevar    resliced
