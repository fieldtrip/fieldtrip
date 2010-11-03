function mri = ft_volumereslice(cfg, mri)

% FT_VOLUMERESLICE reslices a volume along the principal axes of the
% coordinate system according to a specified resolution.
%
% Use as
%   mri = ft_volumereslice(cfg, mri)
% where the mri contains an anatomical or functional volume and cfg is a
% configuration structure containing
%   cfg.xrange     = [min max], in physical units
%   cfg.yrange     = [min max], in physical units
%   cfg.zrange     = [min max], in physical units
%   cfg.resolution = number, in physical units
%
% See also FT_VOLUMEDOWNSAMPLE, FT_SOURCEINTERPOLATE

% Undocumented local options:
%   cfg.inputfile        = one can specifiy preanalysed saved data as input
%   cfg.outputfile       = one can specify output as file to save to disk

% Copyright (C) 2010, Robert Oostenveld
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set the defaults
if ~isfield(cfg, 'resolution');   cfg.resolution   = 1;         end % in physical units
if ~isfield(cfg, 'downsample');   cfg.downsample   = 1;         end
if ~isfield(cfg, 'inputfile'),    cfg.inputfile    = [];        end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile   = [];        end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    mri = loadvar(cfg.inputfile, 'mri');
  end
end

% check if the input data is valid for this function and ensure that the structures correctly describes a volume
mri = ft_checkdata(mri, 'datatype', 'volume', 'inside', 'logical', 'feedback', 'yes', 'hasunits', 'yes');

cfg = ft_checkconfig(cfg, 'required', {'xrange', 'yrange', 'zrange'});

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

pseudomri           = [];
pseudomri.dim       = [length(xgrid) length(ygrid) length(zgrid)];
pseudomri.transform = translate([cfg.xrange(1) cfg.yrange(1) cfg.zrange(1)]) * scale([cfg.resolution cfg.resolution cfg.resolution]);
pseudomri.anatomy   = zeros(pseudomri.dim, 'int8');

clear xgrid ygrid zgrid

fprintf('reslicing from [%d %d %d] to [%d %d %d]\n', mri.dim(1), mri.dim(2), mri.dim(3), pseudomri.dim(1), pseudomri.dim(2), pseudomri.dim(3));

tmpcfg = [];
mri = ft_sourceinterpolate(tmpcfg, mri, pseudomri);

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: ft_sourceinterpolate.m 715 2010-03-09 10:57:27Z roboos $';
% remember the configuration details of the input data
try cfg.previous = mri.cfg; end
% remember the exact configuration details in the output
mri.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'mri', mri); % use the variable name "data" in the output file
end
