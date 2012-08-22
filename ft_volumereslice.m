function [reslice] = ft_volumereslice(cfg, mri)

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
% The new spatial extent can be specified with
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
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
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

% Copyright (C) 2010-2011, Robert Oostenveld & Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar mri

% check if the input data is valid for this function and ensure that the structures correctly describes a volume
mri = ft_checkdata(mri, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunits', 'yes');

% set the defaults
cfg.resolution = ft_getopt(cfg, 'resolution', 1);
cfg.downsample = ft_getopt(cfg, 'downsample', 1);
cfg.inputfile  = ft_getopt(cfg, 'inputfile',  []);
cfg.outputfile = ft_getopt(cfg, 'outputfile', []);
cfg.xrange     = ft_getopt(cfg, 'xrange', []);
cfg.yrange     = ft_getopt(cfg, 'yrange', []);
cfg.zrange     = ft_getopt(cfg, 'zrange', []);

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
      zshift = 15./cfg.resolution;
  end
else
  xshift = 0;
  yshift = 0;
  zshift = 0;
end

cfg.dim = ft_getopt(cfg, 'dim',    ceil(mri.dim./cfg.resolution));
if isempty(cfg.xrange),
  cfg.xrange = [-cfg.dim(1)/2+0.5 cfg.dim(1)/2-0.5] * cfg.resolution + xshift;
end
if isempty(cfg.yrange),
  cfg.yrange = [-cfg.dim(2)/2+0.5 cfg.dim(2)/2-0.5] * cfg.resolution + yshift;
end
if isempty(cfg.zrange),
  cfg.zrange = [-cfg.dim(3)/2+0.5 cfg.dim(3)/2-0.5] * cfg.resolution + zshift;
end

if ~isequal(cfg.downsample, 1)
  % downsample the anatomical volume
  tmpcfg = [];
  tmpcfg.downsample = cfg.downsample;
  mri = ft_volumedownsample(tmpcfg, mri);
end

% compute the desired grid positions
xgrid = cfg.xrange(1):cfg.resolution:cfg.xrange(2);
ygrid = cfg.yrange(1):cfg.resolution:cfg.yrange(2);
zgrid = cfg.zrange(1):cfg.resolution:cfg.zrange(2);

reslice           = [];
reslice.dim       = [length(xgrid) length(ygrid) length(zgrid)];
reslice.transform = translate([cfg.xrange(1) cfg.yrange(1) cfg.zrange(1)]) * scale([cfg.resolution cfg.resolution cfg.resolution]) * translate([-1 -1 -1]);
reslice.anatomy   = zeros(reslice.dim, 'int8');

clear xgrid ygrid zgrid

% these are the same in the resliced as in the input anatomical MRI
if isfield(mri, 'coordsys')
  reslice.coordsys = mri.coordsys;
end
if isfield(mri, 'unit')
  reslice.unit = mri.unit;
end

fprintf('reslicing from [%d %d %d] to [%d %d %d]\n', mri.dim(1), mri.dim(2), mri.dim(3), reslice.dim(1), reslice.dim(2), reslice.dim(3));

% the actual work is being done by ft_sourceinterpolate, which interpolates the real mri volume 
% on the resolution that is defined for the resliced volume
tmpcfg = [];
tmpcfg.parameter = 'anatomy';
reslice = ft_sourceinterpolate(tmpcfg, mri, reslice);

% convert any non-finite values to 0 to avoid problems later on
reslice.anatomy(~isfinite(reslice.anatomy)) = 0;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous mri
ft_postamble history reslice
ft_postamble savevar reslice
