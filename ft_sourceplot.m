function [cfg] = ft_sourceplot(cfg, data)

% FT_SOURCEPLOT plots functional source reconstruction data on slices or on
% a surface, optionally as an overlay on anatomical MRI data, where
% statistical data can be used to determine the opacity of the mask. Input
% data comes from FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE or statistical
% values from FT_SOURCESTATISTICS.
%
% Use as
%   ft_sourceplot(cfg, data)
% where the input data can contain an anatomical MRI, functional source
% reconstruction results and statistical data, which all have to be
% interpolated on the same 3-d grid.
%
% The configuration should contain:
%   cfg.method        = 'slice',   plots the data on a number of slices in the same plane
%                       'ortho',   plots the data on three orthogonal slices
%                       'surface', plots the data on a 3D brain surface
%
%   cfg.anaparameter  = string, field in data with the anatomical data (default = 'anatomy' if present in data)
%   cfg.funparameter  = string, field in data with the functional parameter of interest (default = [])
%   cfg.maskparameter = string, field in the data to be used for opacity masking of fun data (default = [])
%                        If values are between 0 and 1, zero is fully transparant and one is fully opaque.
%                        If values in the field are not between 0 and 1 they will be scaled depending on the values
%                        of cfg.opacitymap and cfg.opacitylim (see below)
%                        You can use masking in several ways, f.i.
%                        - use outcome of statistics to show only the significant values and mask the insignificant
%                          NB see also cfg.opacitymap and cfg.opacitylim below
%                        - use the functional data itself as mask, the highest value (and/or lowest when negative)
%                          will be opaque and the value closest to zero transparent
%                        - Make your own field in the data with values between 0 and 1 to control opacity directly
%
% The following parameters can be used in all methods:
%   cfg.downsample    = downsampling for resolution reduction, integer value (default = 1) (orig: from surface)
%   cfg.atlas         = string, filename of atlas to use (default = []) SEE FT_PREPARE_ATLAS
%                        for ROI masking (see "masking" below) or in interactive mode (see "ortho-plotting" below)
%   cfg.coordsys      = 'mni' or 'tal', coordinate system of the input data, used to lookup the label from the atlas
%
% The following parameters can be used for the functional data:
%   cfg.funcolormap   = colormap for functional data, see COLORMAP (default = 'auto')
%                       'auto', depends structure funparameter, or on funcolorlim
%                         - funparameter: only positive values, or funcolorlim:'zeromax' -> 'hot'
%                         - funparameter: only negative values, or funcolorlim:'minzero' -> 'cool'
%                         - funparameter: both pos and neg values, or funcolorlim:'maxabs' -> 'jet'
%                         - funcolorlim: [min max] if min & max pos-> 'hot', neg-> 'cool', both-> 'jet'
%   cfg.funcolorlim   = color range of the functional data (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(funparameter)) to +max(abs(funparameter))
%                        'zeromax', from 0 to max(abs(funparameter))
%                        'minzero', from min(abs(funparameter)) to 0
%                        'auto', if funparameter values are all positive: 'zeromax',
%                          all negative: 'minzero', both possitive and negative: 'maxabs'
%   cfg.colorbar      = 'yes' or 'no' (default = 'yes')
%
% The following parameters can be used for the masking data:
%   cfg.opacitymap    = opacitymap for mask data, see ALPHAMAP (default = 'auto')
%                       'auto', depends structure maskparameter, or on opacitylim
%                         - maskparameter: only positive values, or opacitylim:'zeromax' -> 'rampup'
%                         - maskparameter: only negative values, or opacitylim:'minzero' -> 'rampdown'
%                         - maskparameter: both pos and neg values, or opacitylim:'maxabs' -> 'vdown'
%                         - opacitylim: [min max] if min & max pos-> 'rampup', neg-> 'rampdown', both-> 'vdown'
%                         - NB. to use p-values use 'rampdown' to get lowest p-values opaque and highest transparent
%   cfg.opacitylim    = range of mask values to which opacitymap is scaled (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(maskparameter)) to +max(abs(maskparameter))
%                        'zeromax', from 0 to max(abs(maskparameter))
%                        'minzero', from min(abs(maskparameter)) to 0
%                        'auto', if maskparameter values are all positive: 'zeromax',
%                          all negative: 'minzero', both possitive and negative: 'maxabs'
%   cfg.roi           = string or cell of strings, region(s) of interest from anatomical atlas (see cfg.atlas above)
%                        everything is masked except for ROI
%
% The folowing parameters apply for ortho-plotting
%   cfg.location      = location of cut, (default = 'auto')
%                        'auto', 'center' if only anatomy, 'max' if functional data
%                        'min' and 'max' position of min/max funparameter
%                        'center' of the brain
%                        [x y z], coordinates in voxels or head, see cfg.locationcoordinates
%   cfg.locationcoordinates = coordinate system used in cfg.location, 'head' or 'voxel' (default = 'head')
%                              'head', headcoordinates from anatomical MRI
%                              'voxel', voxelcoordinates
%   cfg.crosshair     = 'yes' or 'no' (default = 'yes')
%   cfg.axis          = 'on' or 'off' (default = 'on')
%   cfg.interactive   = 'yes' or 'no' (default = 'no')
%                        in interactive mode cursor click determines location of cut
%   cfg.queryrange    = number, in atlas voxels (default 3)
%
%
% The folowing parameters apply for slice-plotting
%   cfg.nslices       = number of slices, (default = 20)
%   cfg.slicerange    = range of slices in data, (default = 'auto')
%                       'auto', full range of data
%                       [min max], coordinates of first and last slice in voxels
%   cfg.slicedim      = dimension to slice 1 (x-axis) 2(y-axis) 3(z-axis) (default = 3)
%   cfg.title         = string, title of the figure window
%
% The folowing parameters apply for surface-plotting
%   cfg.surffile       = string, file that contains the surface (default = 'single_subj_T1.mat')
%                        'single_subj_T1.mat' contains a triangulation that corresponds with the
%                         SPM anatomical template in MNI coordinates
%   cfg.surfinflated   = string, file that contains the inflated surface (default = [])
%   cfg.surfdownsample = number (default = 1, i.e. no downsampling)
%   cfg.projmethod     = projection method, how functional volume data is
%   projected onto surface
%                        'nearest', 'project', 'sphere_avg', 'sphere_weighteddistance'
%   cfg.projvec        = vector (in mm) to allow different projections that
%   are combined with the method specified in cfg.projcomb
%   cfg.projcomb       = 'mean', 'max', method to combine the different
% projections
%   cfg.projweight     = vector of weights for the different projections
%   (default: 1)
%   cfg.projthresh      = implements thresholding on the surface level
%   (cfg.projthresh = 0.7 means 70% of maximum)
%   cfg.sphereradius   = maximum distance from each voxel to the surface to be
%                        included in the sphere projection methods, expressed in mm
%   cfg.distmat        = precomputed distance matrix (default = [])
%   cfg.camlight       = 'yes' or 'no' (default = 'yes')
%   cfg.renderer       = 'painters', 'zbuffer',' opengl' or 'none' (default = 'opengl')
%                        When using opacity the OpenGL renderer is required.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following option:
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE, FT_SOURCESTATISTICS,
%  FT_VOLUMELOOKUP, FT_PREPARE_ATLAS

% TODO have to be built in:
%   cfg.marker        = [Nx3] array defining N marker positions to display (orig: from sliceinterp)
%   cfg.markersize    = radius of markers (default = 5)
%   cfg.markercolor   = [1x3] marker color in RGB (default = [1 1 1], i.e. white) (orig: from sliceinterp)
%   white background option

% undocumented TODO
%   slice in all directions
%   surface also optimal when inside present
%   come up with a good glass brain projection

% Copyright (C) 2007-2008, Robert Oostenveld, Ingrid Nieuwenhuis
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
ft_preamble loadvar data

% this is not supported any more as of 26/10/2011
if ischar(data)
  error('please use cfg.inputfile instead of specifying the input variable as a sting');
end

% check if the input data is valid for this function
data     = ft_checkdata(data, 'datatype', {'volume' 'source'}, 'feedback', 'yes');
issource = ft_datatype(data, 'source');
isvolume = ft_datatype(data, 'volume');
if issource && ~isfield(data, 'dim')
  error('the input data needs to be defined on a regular 3D grid');
end

% check if the input configuration is valid for this function
cfg = ft_checkconfig(cfg, 'renamed', {'inputcoordsys', 'coordsys'});

% set the defaults for all methods
cfg.method        = ft_getopt(cfg, 'method', 'ortho');
cfg.funparameter  = ft_getopt(cfg, 'funparameter',  []);
cfg.maskparameter = ft_getopt(cfg, 'maskparameter', []);
cfg.downsample    = ft_getopt(cfg, 'downsample',    1);
cfg.title         = ft_getopt(cfg, 'title',         '');
cfg.atlas         = ft_getopt(cfg, 'atlas',         []);
cfg.marker        = ft_getopt(cfg, 'marker',        []);
cfg.markersize    = ft_getopt(cfg, 'markersize',    5);
cfg.markercolor   = ft_getopt(cfg, 'markercolor',   [1 1 1]);

if ~isfield(cfg, 'anaparameter')
  if isfield(data, 'anatomy')
    cfg.anaparameter = 'anatomy';
  else
    cfg.anaparameter = [];
  end
end

% set the common defaults for the functional data
cfg.funcolormap   = ft_getopt(cfg, 'funcolormap',   'auto');
cfg.funcolorlim   = ft_getopt(cfg, 'funcolorlim',   'auto');

% set the common defaults for the statistical data
cfg.opacitymap    = ft_getopt(cfg, 'opacitymap',    'auto');
cfg.opacitylim    = ft_getopt(cfg, 'opacitylim',    'auto');
cfg.roi           = ft_getopt(cfg, 'roi',           []);

% set the defaults per method

% ortho
cfg.location            = ft_getopt(cfg, 'location',            'auto');
cfg.locationcoordinates = ft_getopt(cfg, 'locationcoordinates', 'head');
cfg.crosshair           = ft_getopt(cfg, 'crosshair',           'yes'); 
cfg.colorbar            = ft_getopt(cfg, 'colorbar',            'yes'); 
cfg.axis                = ft_getopt(cfg, 'axis',                'on');  
cfg.interactive         = ft_getopt(cfg, 'interactive',         'no');  
cfg.queryrange          = ft_getopt(cfg, 'queryrange',          3);     
cfg.coordsys            = ft_getopt(cfg, 'coordsys',            []);    
cfg.units               = ft_getopt(cfg, 'units',               []);    

if isfield(cfg, 'TTlookup'),
  error('TTlookup is old; now specify cfg.atlas, see help!');
end

% slice
cfg.nslices    = ft_getopt(cfg, 'nslices',    20);
cfg.slicedim   = ft_getopt(cfg, 'slicedim',   3);
cfg.slicerange = ft_getopt(cfg, 'slicerange', 'auto');

% surface
cfg.downsample     = ft_getopt(cfg, 'downsample',     1);
cfg.surfdownsample = ft_getopt(cfg, 'surfdownsample', 1);
cfg.surffile       = ft_getopt(cfg, 'surffile',       'single_subj_T1.mat');% use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
cfg.surfinflated   = ft_getopt(cfg, 'surfinflated',  []);
cfg.sphereradius   = ft_getopt(cfg, 'sphereradius',  []);
cfg.projvec        = ft_getopt(cfg, 'projvec',       1); 
cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); %or max
cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);                 
cfg.distmat        = ft_getopt(cfg, 'distmat',       []);
cfg.camlight       = ft_getopt(cfg, 'camlight',      'yes');
cfg.renderer       = ft_getopt(cfg, 'renderer',      'opengl');
if isequal(cfg.method,'surface')
  if ~isfield(cfg, 'projmethod'), error('specify cfg.projmethod'); end 
end

% for backward compatibility
if strcmp(cfg.location, 'interactive')
  cfg.location = 'auto';
  cfg.interactive = 'yes';
end

% ensure that the data has interpretable spatial units
if     ~isfield(data, 'unit') && ~isempty(cfg.units)
  data.unit = cfg.units;
elseif ~isfield(data, 'unit') &&  isempty(cfg.units)
  data = ft_convert_units(data);
elseif  isfield(data, 'unit') && ~isempty(cfg.units)
  data = ft_convert_units(data, cfg.units);
elseif  isfield(data, 'unit') &&  isempty(cfg.units)
  % nothing to do
end

% ensure that the data has an interpretable coordinate system
if     ~isfield(data, 'coordsys') && ~isempty(cfg.coordsys)
  data.coordsys = cfg.coordsys;
elseif ~isfield(data, 'coordsys') &&  isempty(cfg.coordsys) && ~isempty(cfg.atlas)
  % only needed if an atlas was specified for volumelookup
  data = ft_convert_coordsys(data);
elseif  isfield(data, 'coordsys') && ~isempty(cfg.coordsys) && ~isempty(cfg.atlas)
  % only needed if an atlas was specified for volumelookup
  data = ft_convert_coordsys(data, cfg.units);
elseif  isfield(data, 'coordsys') &&  isempty(cfg.coordsys)
  % nothing to do
end

% select the functional and the mask parameter
cfg.funparameter  = parameterselection(cfg.funparameter, data);
cfg.maskparameter = parameterselection(cfg.maskparameter, data);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1};  end
try, cfg.maskparameter = cfg.maskparameter{1}; end

if cfg.downsample ~=1 && isvolume
  % downsample all volumes
  tmpcfg = [];
  tmpcfg.parameter  = {cfg.funparameter, cfg.maskparameter, cfg.anaparameter};
  tmpcfg.downsample = cfg.downsample;
  data = ft_volumedownsample(tmpcfg, data);
end

%%% make the local variables:
dim = data.dim;

hasatlas = ~isempty(cfg.atlas);
if hasatlas
  % initialize the atlas
  [p, f, x] = fileparts(cfg.atlas);
  fprintf(['reading ', f,' atlas coordinates and labels\n']);
  atlas = ft_prepare_atlas(cfg);
end

hasroi = ~isempty(cfg.roi);
if hasroi
  if ~hasatlas
    error('specify cfg.atlas which belongs to cfg.roi')
  else
    % get the mask
    tmpcfg          = [];
    tmpcfg.roi      = cfg.roi;
    tmpcfg.atlas    = cfg.atlas;
    tmpcfg.coordsys = cfg.coordsys;
    roi = ft_volumelookup(tmpcfg,data);
  end
end

%%% anaparameter
if isempty(cfg.anaparameter);
  hasana = 0;
  fprintf('not plotting anatomy\n');
elseif isfield(data, cfg.anaparameter)
  hasana = 1;
  ana = getsubfield(data, cfg.anaparameter);
  % convert integers to single precision float if neccessary
  if isa(ana, 'uint8') || isa(ana, 'uint16') || isa(ana, 'int8') || isa(ana, 'int16')
    fprintf('converting anatomy to double\n');
    ana = double(ana);
  end
  fprintf('scaling anatomy to [0 1]\n');
  dmin = min(ana(:));
  dmax = max(ana(:));
  ana  = (ana-dmin)./(dmax-dmin);
else
  warning('do not understand cfg.anaparameter, not plotting anatomy\n')
  hasana = 0;
end

%%% funparameter
% has fun?
if ~isempty(cfg.funparameter)
  if issubfield(data, cfg.funparameter)
    hasfun = 1;
    fun = getsubfield(data, cfg.funparameter);
  else
    error('cfg.funparameter not found in data');
  end
else
  hasfun = 0;
  fprintf('no functional parameter\n');
end

% handle fun
if hasfun && issubfield(data, 'dimord') && strcmp(data.dimord(end-2:end),'rgb')
  % treat functional data as rgb values
  if any(fun(:)>1 | fun(:)<0)
    %scale
    tmpdim = size(fun);
    nvox   = prod(tmpdim(1:end-1));
    tmpfun = reshape(fun,[nvox tmpdim(end)]);
    m1     = max(tmpfun,[],1);
    m2     = min(tmpfun,[],1);
    tmpfun = (tmpfun-m2(ones(nvox,1),:))./(m1(ones(nvox,1),:)-m2(ones(nvox,1),:));
    fun = reshape(tmpfun, tmpdim);
  end
  qi      = 1;
  hasfreq = 0;
  hastime = 0;
  
  doimage = 1;
  fcolmin = 0;
  fcolmax = 1;
elseif hasfun
  % determine scaling min and max (fcolmin fcolmax) and funcolormap
  if ~isa(fun, 'logical')
    funmin = min(fun(:));
    funmax = max(fun(:));
  else
    funmin = 0;
    funmax = 1;
  end
  % smart lims: make from auto other string
  if isequal(cfg.funcolorlim,'auto')
    if sign(funmin)>-1 && sign(funmax)>-1
      cfg.funcolorlim = 'zeromax';
    elseif sign(funmin)<1 && sign(funmax)<1
      cfg.funcolorlim = 'minzero';
    else
      cfg.funcolorlim = 'maxabs';
    end
  end
  if ischar(cfg.funcolorlim)
    % limits are given as string
    if isequal(cfg.funcolorlim,'maxabs')
      fcolmin = -max(abs([funmin,funmax]));
      fcolmax =  max(abs([funmin,funmax]));
      if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'jet'; end;
    elseif isequal(cfg.funcolorlim,'zeromax')
      fcolmin = 0;
      fcolmax = funmax;
      if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'hot'; end;
    elseif isequal(cfg.funcolorlim,'minzero')
      fcolmin = funmin;
      fcolmax = 0;
      if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'cool'; end;
    else
      error('do not understand cfg.funcolorlim');
    end
  else
    % limits are numeric
    fcolmin = cfg.funcolorlim(1);
    fcolmax = cfg.funcolorlim(2);
    % smart colormap
    if isequal(cfg.funcolormap,'auto')
      if sign(fcolmin) == -1 && sign(fcolmax) == 1
        cfg.funcolormap = 'jet';
      else
        if fcolmin < 0
          cfg.funcolormap = 'cool';
        else
          cfg.funcolormap = 'hot';
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

  %what if fun is 4D?
  if ndims(fun)>3 || prod(dim)==size(fun,1)
    if isfield(data, 'time') && isfield(data, 'freq'),
      %data contains timefrequency representation
      qi      = [1 1];
      hasfreq = 1;
      hastime = 1;
      fun     = reshape(fun, [dim numel(data.freq) numel(data.time)]);
    elseif isfield(data, 'time')
      %data contains evoked field
      qi      = 1;
      hasfreq = 0;
      hastime = 1;
      fun     = reshape(fun, [dim numel(data.time)]);
    elseif isfield(data, 'freq')
      %data contains frequency spectra
      qi      = 1;
      hasfreq = 1;
      hastime = 0;
      fun     = reshape(fun, [dim numel(data.freq)]);
    else
      qi      = 1;
      hasfreq = 0;
      hastime = 0;
      fun     = reshape(fun, dim);
    end
  else
    %do nothing
    qi      = 1;
    hasfreq = 0;
    hastime = 0;
  end
  
  doimage = 0;
else
  qi      = 1;
  hasfreq = 0;
  hastime = 0;

  doimage = 0;
end % handle fun

%%% maskparameter
% has mask?
if ~isempty(cfg.maskparameter)
  if issubfield(data, cfg.maskparameter)
    if ~hasfun
      error('you can not have a mask without functional data')
    else
      hasmsk = 1;
      msk = getsubfield(data, cfg.maskparameter);
      if islogical(msk) %otherwise sign() not posible
        msk = double(msk);
      end
    end
  else
    error('cfg.maskparameter not found in data');
  end
else
  hasmsk = 0;
  fprintf('no masking parameter\n');
end
% handle mask
if hasmsk
  % determine scaling and opacitymap
  mskmin = min(msk(:));
  mskmax = max(msk(:));
  % determine the opacity limits and the opacity map
  % smart lims: make from auto other string, or equal to funcolorlim if funparameter == maskparameter
  if isequal(cfg.opacitylim,'auto')
    if isequal(cfg.funparameter,cfg.maskparameter)
      cfg.opacitylim = cfg.funcolorlim;
    else
      if sign(mskmin)>-1 && sign(mskmax)>-1
        cfg.opacitylim = 'zeromax';
      elseif sign(mskmin)<1 && sign(mskmax)<1
        cfg.opacitylim = 'minzero';
      else
        cfg.opacitylim = 'maxabs';
      end
    end
  end
  if ischar(cfg.opacitylim)
    % limits are given as string
    switch cfg.opacitylim
      case 'zeromax'
        opacmin = 0;
        opacmax = mskmax;
        if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'rampup'; end;
      case 'minzero'
        opacmin = mskmin;
        opacmax = 0;
        if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'rampdown'; end;
      case 'maxabs'
        opacmin = -max(abs([mskmin, mskmax]));
        opacmax =  max(abs([mskmin, mskmax]));
        if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'vdown'; end;
      otherwise
        error('incorrect specification of cfg.opacitylim');
    end % switch opacitylim
  else
    % limits are numeric
    opacmin = cfg.opacitylim(1);
    opacmax = cfg.opacitylim(2);
    if isequal(cfg.opacitymap,'auto')
      if sign(opacmin)>-1 && sign(opacmax)>-1
        cfg.opacitymap = 'rampup';
      elseif sign(opacmin)<1 && sign(opacmax)<1
        cfg.opacitymap = 'rampdown';
      else
        cfg.opacitymap = 'vdown';
      end
    end
  end % handling opacitylim and opacitymap
  clear mskmin mskmax;
end

% prevent outside fun from being plotted
if hasfun && isfield(data,'inside') && ~hasmsk
  hasmsk = 1;
  msk = zeros(dim);
  cfg.opacitymap = 'rampup';
  opacmin = 0;
  opacmax = 1;
  % make intelligent mask
  if isequal(cfg.method,'surface')
    msk(data.inside) = 1;
  else
    if hasana
      msk(data.inside) = 0.5; %so anatomy is visible
    else
      msk(data.inside) = 1;
    end
  end
end

% if region of interest is specified, mask everything besides roi
if hasfun && hasroi && ~hasmsk
  hasmsk = 1;
  msk = roi;
  cfg.opacitymap = 'rampup';
  opacmin = 0;
  opacmax = 1;
elseif hasfun && hasroi && hasmsk
  msk = roi .* msk;
elseif hasroi
  error('you can not have a roi without functional data')
end

%%% set color and opacity mapping for this figure
if hasfun
  cfg.funcolormap = colormap(cfg.funcolormap);
  colormap(cfg.funcolormap);
end
if hasmsk
  cfg.opacitymap  = alphamap(cfg.opacitymap);
  alphamap(cfg.opacitymap);
  if ndims(fun)>3 && ndims(msk)==3
    siz = size(fun);
    msk = repmat(msk, [1 1 1 siz(4:end)]);
  end
end

%%% determine what has to be plotted, depends on method
if isequal(cfg.method,'ortho')
  if ~ischar(cfg.location)
    if strcmp(cfg.locationcoordinates, 'head')
      % convert the headcoordinates location into voxel coordinates
      loc = inv(data.transform) * [cfg.location(:); 1];
      loc = round(loc(1:3));
    elseif strcmp(cfg.locationcoordinates, 'voxel')
      % the location is already in voxel coordinates
      loc = round(cfg.location(1:3));
    else
      error('you should specify cfg.locationcoordinates');
    end
  else
    if isequal(cfg.location,'auto')
      if hasfun
        if isequal(cfg.funcolorlim,'maxabs');
          loc = 'max';
        elseif isequal(cfg.funcolorlim, 'zeromax');
          loc = 'max';
        elseif isequal(cfg.funcolorlim, 'minzero');
          loc = 'min';
        else %if numerical
          loc = 'max';
        end
      else
        loc = 'center';
      end;
    else
      loc = cfg.location;
    end
  end

  % determine the initial intersection of the cursor (xi yi zi)
  if ischar(loc) && strcmp(loc, 'min')
    if isempty(cfg.funparameter)
      error('cfg.location is min, but no functional parameter specified');
    end
    [~, minindx] = min(fun(:));
    [xi, yi, zi] = ind2sub(dim, minindx);
  elseif ischar(loc) && strcmp(loc, 'max')
    if isempty(cfg.funparameter)
      error('cfg.location is max, but no functional parameter specified');
    end
    [~, maxindx] = max(fun(:));
    [xi, yi, zi] = ind2sub(dim, maxindx);
  elseif ischar(loc) && strcmp(loc, 'center')
    xi = round(dim(1)/2);
    yi = round(dim(2)/2);
    zi = round(dim(3)/2);
  elseif ~ischar(loc)
    % using nearest instead of round ensures that the position remains within the volume
    xi = nearest(1:dim(1), loc(1));
    yi = nearest(1:dim(2), loc(2));
    zi = nearest(1:dim(3), loc(3));
  end

  % enforce the size of the subplots to be isotropic
  xdim = dim(1) + dim(2);
  ydim = dim(2) + dim(3);
  
  xsize(1) = 0.82*dim(1)/xdim;
  xsize(2) = 0.82*dim(2)/xdim;
  ysize(1) = 0.82*dim(3)/ydim;
  ysize(2) = 0.82*dim(2)/ydim;
 
  %% do the actual plotting %%
  nas = [];
  lpa = [];
  rpa = [];
  tag = 'ik';
  update = [1 1 1];
  interactive_flag = 1; % it happens at least once
  feedbackmessage  = true;
  while(interactive_flag)
    interactive_flag = strcmp(cfg.interactive, 'yes');

    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

    if interactive_flag && feedbackmessage
      fprintf('\n');
      fprintf('click with mouse button to reposition the cursor\n');
      fprintf('press n/l/r on keyboard to record a fiducial position\n');
      fprintf('press q on keyboard to quit interactive mode\n');
      feedbackmessage = false;
    end
    
    ijk = [xi yi zi 1]';
    if isvolume
      xyz = data.transform * ijk;
    elseif issource
      ix  = sub2ind(dim,xi,yi,zi);
      xyz = data.pos(ix,:);
    end
    ijk = ijk(1:3);
    
    % construct a string with user feedback
    str1 = sprintf('voxel %d, indices [%d %d %d]', sub2ind(dim(1:3), xi, yi, zi), ijk);

    if isfield(data, 'coordsys') && isfield(data, 'unit')
      str2 = sprintf('%s coordinates [%.1f %.1f %.1f] %s', data.coordsys, xyz(1:3), data.unit);
    elseif ~isfield(data, 'coordsys') && isfield(data, 'unit')
      str2 = sprintf('location [%.1f %.1f %.1f] %s', xyz(1:3), data.unit);
    elseif isfield(data, 'coordsys') && ~isfield(data, 'unit')
      str2 = sprintf('%s coordinates [%.1f %.1f %.1f]', data.coordsys, xyz(1:3));
    elseif ~isfield(data, 'coordsys') && ~isfield(data, 'unit')
      str2 = sprintf('location [%.1f %.1f %.1f]', xyz(1:3));
    else
      str2 = '';
    end
    
    if hasfreq && hastime,
      str3 = sprintf('%.1f s, %.1f Hz', data.freq(qi(1)), data.time(qi(2)));
    elseif ~hasfreq && hastime,
      str3 = sprintf('%.1f s', data.time(qi(1)));
    elseif hasfreq && ~hastime,
      str3 = sprintf('%.1f Hz', data.freq(qi(1)));
    else 
      str3 = '';
    end
    
    if hasfun
      if ~hasfreq && ~hastime
        val = fun(xi, yi, zi);
      elseif ~hasfreq && hastime
        val = fun(xi, yi, zi, qi);
      elseif hasfreq && ~hastime
        val = fun(xi, yi, zi, qi);
      elseif hasfreq && hastime
        val = fun(xi, yi, zi, qi(1), qi(2));
      end
      str4 = sprintf('value %f', val);
    else
      str4 = '';
    end
    
    %fprintf('%s %s %s %s\n', str1, str2, str3, str4);

    if hasatlas
      % determine the anatomical label of the current position
      lab = atlas_lookup(atlas, (xyz(1:3)), 'inputcoord', data.coordsys, 'queryrange', cfg.queryrange);
      if isempty(lab)
        fprintf([f,' labels: not found\n']);
      else
        fprintf([f,' labels: '])
        fprintf('%s', lab{1});
        for i=2:length(lab)
          fprintf(', %s', lab{i});
        end
        fprintf('\n');
      end
    end
    
    if ~hasfun && ~hasana
      % this seems to be a problem that people often have
      error('no anatomy is present and no functional data is selected, please check your cfg.funparameter');
    end
    
    % create only once
    if ~exist('h1a','var')
      % these will hold the anatomical data if present, along with labels
      % etc.
      h1a = axes('position',[0.07 0.07+ysize(2)+0.05 xsize(1) ysize(1)]);
      h2a = axes('position',[0.07+xsize(1)+0.05 0.07+ysize(2)+0.05 xsize(2) ysize(1)]);
      h3a = axes('position',[0.07 0.07 xsize(1) ysize(2)]);
      set(h1a,'Tag','ik','Visible',cfg.axis,'XAxisLocation','top');
      set(h2a,'Tag','jk','Visible',cfg.axis,'XAxisLocation','top');
      set(h3a,'Tag','ij','Visible',cfg.axis);
      set(gcf, 'renderer', cfg.renderer); % ensure that this is done in interactive mode
    end
      
    if hasana,
      ft_plot_ortho(ana, 'transform', eye(4), 'location', ijk, 'style', 'subplot', 'parents', [h1a h2a h3a].*update, 'doscale', false);
    end
    if hasfun,
      % create only once: these are different axes than the anatomy -> could be different resolution (to be implemented)
      if ~exist('h1f','var')
        h1f = axes('position',[0.07 0.07+ysize(2)+0.05 xsize(1) ysize(1)]);
        h2f = axes('position',[0.07+xsize(1)+0.05 0.07+ysize(2)+0.05 xsize(2) ysize(1)]);
        h3f = axes('position',[0.07 0.07 xsize(1) ysize(2)]);
        set(h1f,'Tag','ik','visible','off');
        set(h2f,'Tag','jk','visible','off');
        set(h3f,'Tag','ij','visible','off');
      end
      if hasmsk
        ft_plot_ortho(fun(:,:,:,qi), msk(:,:,:,qi), 'transform', eye(4), 'location', ijk, ...
                           'style', 'subplot', 'parents', [h1f h2f h3f].*update, ...
                           'colormap', cfg.funcolormap, 'colorlim', [fcolmin fcolmax], ...
                           'opacitylim', [opacmin opacmax]);
      else
        ft_plot_ortho(fun(:,:,:,qi), 'transform', eye(4), 'location', ijk, ...
                           'style', 'subplot', 'parents', [h1f h2f h3f].*update, ...
                           'colormap', cfg.funcolormap, 'colorlim', [fcolmin fcolmax]);
      end
    end
   
    maketransparent = true;
    if maketransparent && ~hasmsk && hasfun && hasana
      set(findobj(get(h1f,'children'),'type','surface'),'facealpha',0.5);
      set(findobj(get(h2f,'children'),'type','surface'),'facealpha',0.5);
      set(findobj(get(h3f,'children'),'type','surface'),'facealpha',0.5);
      maketransparent = false;
    end
    
%     if strcmp(cfg.crosshair, 'yes'),
%       if update(1), set(gcf,'currentaxes',h1f); crosshair([xi zi]); end
%       if update(2), set(gcf,'currentaxes',h2f); crosshair([yi zi]); end
%       if update(3), set(gcf,'currentaxes',h3f); crosshair([xi yi]); end
%     end
%     
       
    if hasfreq && hastime && hasfun,
      h=subplot(2,2,4);
      %uimagesc(data.time, data.freq, squeeze(vols{2}(xi,yi,zi,:,:))');axis xy;
      tmpdat = double(squeeze(fun(xi,yi,zi,:,:)));
      pcolor(double(data.time), double(data.freq), double(squeeze(vols{2}(xi,yi,zi,:,:))));
      shading('interp');
      xlabel('time'); ylabel('freq');
      try
        caxis([-1 1].*max(abs(caxis)));
      end
      set(gca,'tag','TF1');
      colorbar;
      %caxis([fcolmin fcolmax]);
      %set(gca, 'Visible', 'off');
    elseif hasfreq && hasfun,
      subplot(2,2,4);
      plot(data.freq, squeeze(fun(xi,yi,zi,:))); xlabel('freq');
      axis([data.freq(1) data.freq(end) fcolmin fcolmax]);
      set(gca,'tag','TF2');
    elseif hastime && hasfun,
      subplot(2,2,4);
      plot(data.time, squeeze(fun(xi,yi,zi,:))); xlabel('time');
      axis([data.time(1) data.time(end) fcolmin fcolmax]);
      set(gca,'tag','TF3');
    elseif strcmp(cfg.colorbar,  'yes') && ~exist('hc', 'var'),
      if hasfun
        % vectorcolorbar = linspace(fcolmin, fcolmax,length(cfg.funcolormap));
        % imagesc(vectorcolorbar,1,vectorcolorbar);colormap(cfg.funcolormap);
        % use a normal Matlab colorbar, attach it to the invisible 4th subplot
        try
          caxis([fcolmin fcolmax]);
        end
        hc = colorbar;
        set(hc,'location','southoutside');
        set(hc,'position',[0.07+xsize(1)+0.05 0.07+ysize(2)-0.05 xsize(2) 0.05]);
        
        try
          set(hc, 'XLim', [fcolmin fcolmax]);
        end
      else
        warning_once('no colorbar possible without functional data');
      end
    end

    if ~(hasfreq || hastime) && ~exist('ht1')%,'var')
      subplot('position',[0.07+xsize(1)+0.05 0.07 xsize(2) ysize(2)]);
      set(gca,'visible','off');
      ht1=text(0,0.6,str1);
      ht2=text(0,0.5,str2);
      ht3=text(0,0.4,str4);
      ht4=text(0,0.3,str3);
    elseif ~(hasfreq || hastime)
      set(ht1,'string',str1);
      set(ht2,'string',str2);
      set(ht3,'string',str4);
      set(ht4,'string',str3);
    end
      
    % make the last current axes current again
    sel = findobj('type','axes','tag',tag);
    if ~isempty(sel)
      set(gcf, 'currentaxes', sel(1));
    end
      
    if interactive_flag

      try
        kk = waitforbuttonpress;
      catch
        % this happens if the figure is closed
        key='q';
      end

      if kk==0
        % mouse press
        key = '';
      else
        % key press get the value for the key
        key = get(gcf, 'currentcharacter');
      end
      
      tag = get(gca, 'Tag');
      if ~isempty(tag) && kk==0
        ijk = mean(get(gca,'currentpoint'));
        if strcmp(tag, 'ik')
          xi  = round(ijk(1));
          zi  = round(ijk(3));
        elseif strcmp(tag, 'ij')
          xi  = round(ijk(1));
          yi  = round(ijk(2));
        elseif strcmp(tag, 'jk')
          yi  = round(ijk(2));
          zi  = round(ijk(3));
        elseif strcmp(tag, 'TF1')
          % timefreq 
          qi(1) = nearest(data.time, ijk(1));
          qi(2) = nearest(data.freq, ijk(2));
        elseif strcmp(tag, 'TF2')
          % freq only
          qi  = nearest(data.freq, ijk(1));
        elseif strcmp(tag, 'TF3')
          % time only
          qi  = nearest(data.time, ijk(1));
        end
        update = [1 1 1];
      end
      
      if isempty(key)
        % this happens if you press the apple key
        key = '';
      end
      switch key
        case ''
          % do nothing
        case 'q'
          break;
        case 'l'
          lpa = [xi yi zi];
        case 'r'
          rpa = [xi yi zi];
        case 'n'
          nas = [xi yi zi];
        case {'i' 'j' 'k' 'm' 28 29 30 31}
          % update the view to a new position
          if     strcmp(tag,'ik') && (strcmp(key,'i') || key==30), zi = zi+1; update = [0 0 1];
          elseif strcmp(tag,'ik') && (strcmp(key,'j') || key==28), xi = xi-1; update = [0 1 0];
          elseif strcmp(tag,'ik') && (strcmp(key,'k') || key==29), xi = xi+1; update = [0 1 0];
          elseif strcmp(tag,'ik') && (strcmp(key,'m') || key==31), zi = zi-1; update = [0 0 1];
          elseif strcmp(tag,'ij') && (strcmp(key,'i') || key==30), yi = yi+1; update = [1 0 0];
          elseif strcmp(tag,'ij') && (strcmp(key,'j') || key==28), xi = xi-1; update = [0 1 0];
          elseif strcmp(tag,'ij') && (strcmp(key,'k') || key==29), xi = xi+1; update = [0 1 0];
          elseif strcmp(tag,'ij') && (strcmp(key,'m') || key==31), yi = yi-1; update = [1 0 0];
          elseif strcmp(tag,'jk') && (strcmp(key,'i') || key==30), zi = zi+1; update = [0 1 0];
          elseif strcmp(tag,'jk') && (strcmp(key,'j') || key==28), yi = yi-1; update = [1 0 0];
          elseif strcmp(tag,'jk') && (strcmp(key,'k') || key==29), yi = yi+1; update = [1 0 0];
          elseif strcmp(tag,'jk') && (strcmp(key,'m') || key==31), zi = zi-1; update = [1 0 0];
          else
            % do nothing
          end;
        otherwise
        
      end % switch key
    end % if interactive_flag

    
    %if ~isempty(nas), fprintf('nas = [%f %f %f]\n', nas); cfg.fiducial.nas = nas; else fprintf('nas = undefined\n'); end
    %if ~isempty(lpa), fprintf('lpa = [%f %f %f]\n', lpa); cfg.fiducial.lpa = lpa; else fprintf('lpa = undefined\n'); end
    %if ~isempty(rpa), fprintf('rpa = [%f %f %f]\n', rpa); cfg.fiducial.rpa = rpa; else fprintf('rpa = undefined\n'); end
  end % while interactive_flag

elseif isequal(cfg.method,'glassbrain')
  tmpcfg          = [];
  tmpcfg.funparameter = cfg.funparameter;
  tmpcfg.method   = 'ortho';
  tmpcfg.location = [1 1 1];
  tmpcfg.funcolorlim = cfg.funcolorlim;
  tmpcfg.funcolormap = cfg.funcolormap;
  tmpcfg.opacitylim  = cfg.opacitylim;
  tmpcfg.locationcoordinates = 'voxel';
  tmpcfg.maskparameter       = 'inside';
  tmpcfg.axis                = cfg.axis;
  tmpcfg.renderer            = cfg.renderer;
  if hasfun,
    fun = getsubfield(data, cfg.funparameter);
    fun(1,:,:) = max(fun, [], 1);
    fun(:,1,:) = max(fun, [], 2);
    fun(:,:,1) = max(fun, [], 3);
    data = setsubfield(data, cfg.funparameter, fun);
  end

  if hasana,
    ana = getsubfield(data, cfg.anaparameter);
    %ana(1,:,:) = max(ana, [], 1);
    %ana(:,1,:) = max(ana, [], 2);
    %ana(:,:,1) = max(ana, [], 3);
    data = setsubfield(data, cfg.anaparameter, ana);
  end

  if hasmsk,
    msk = getsubfield(data, 'inside');
    msk(1,:,:) = squeeze(fun(1,:,:))>0 & imfill(abs(squeeze(ana(1,:,:))-1))>0;
    msk(:,1,:) = squeeze(fun(:,1,:))>0 & imfill(abs(squeeze(ana(:,1,:))-1))>0;
    msk(:,:,1) = squeeze(fun(:,:,1))>0 & imfill(abs(ana(:,:,1)-1))>0;
    data = setsubfield(data, 'inside', msk);
  end

  ft_sourceplot(tmpcfg, data);

elseif isequal(cfg.method,'surface')

  % read the triangulated cortical surface from file
  tmp = load(cfg.surffile, 'bnd');
  surf = tmp.bnd;
  if isfield(surf, 'transform'),
    % compute the surface vertices in head coordinates
    surf.pnt = warp_apply(surf.transform, surf.pnt);
  end

  % downsample the cortical surface
  if cfg.surfdownsample > 1
    if ~isempty(cfg.surfinflated)
      error('downsampling the surface is not possible in combination with an inflated surface');
    end
    fprintf('downsampling surface from %d vertices\n', size(surf.pnt,1));
    [surf.tri, surf.pnt] = reducepatch(surf.tri, surf.pnt, 1/cfg.surfdownsample);
  end

  % these are required
  if ~isfield(data, 'inside')
    data.inside = true(dim);
  end

  fprintf('%d voxels in functional data\n', prod(dim));
  fprintf('%d vertices in cortical surface\n', size(surf.pnt,1));
  
  if (hasfun  && strcmp(cfg.projmethod,'project')),
	  val=zeros(size(surf.pnt,1),1);
	  if hasmsk
		  maskval = val;
	  end;
	  %convert projvec in mm to a factor, assume mean distance of 70mm
	  cfg.projvec=(70-cfg.projvec)/70;
	  for iproj = 1:length(cfg.projvec),
		  sub = round(warp_apply(inv(data.transform), surf.pnt*cfg.projvec(iproj), 'homogenous'));  % express
		  sub(sub(:)<1) = 1;
		  sub(sub(:,1)>dim(1),1) = dim(1);
		  sub(sub(:,2)>dim(2),2) = dim(2);
		  sub(sub(:,3)>dim(3),3) = dim(3);
		  disp('projecting...')
		  ind = sub2ind(dim, sub(:,1), sub(:,2), sub(:,3));
		  if strcmp(cfg.projcomb,'mean')
			  val = val + cfg.projweight(iproj) * fun(ind);
			  if hasmsk
				  maskval = maskval + cfg.projweight(iproj) * msk(ind);
			  end
		  elseif strcmp(cfg.projcomb,'max')
			  val =  max([val cfg.projweight(iproj) * fun(ind)],[],2);
			  tmp2 = min([val cfg.projweight(iproj) * fun(ind)],[],2);
			  fi = find(val < max(tmp2));
			  val(fi) = tmp2(fi);
			  if hasmsk
				  maskval = max(abs([maskval cfg.projweight(iproj) * fun(ind)]),[],2);
			  end
		  else
			  error('undefined method to combine projections; use cfg.projcomb= mean or max')
		  end
	  end
	  if strcmp(cfg.projcomb,'mean'),
		  val=val/length(cfg.projvec);
		  if hasmsk
			  maskval = max(abs([maskval cfg.projweight(iproj) * fun(ind)]),[],2);
		  end
	  end;
	  if ~isempty(cfg.projthresh),
		  mm=max(abs(val(:)));
		  maskval(abs(val) < cfg.projthresh*mm) = 0;
	  end
  end
  
  if (hasfun && ~strcmp(cfg.projmethod,'project')),
	  [interpmat, cfg.distmat] = interp_gridded(data.transform, fun, surf.pnt, 'projmethod', cfg.projmethod, 'distmat', cfg.distmat, 'sphereradius', cfg.sphereradius, 'inside', data.inside);
	  % interpolate the functional data
	  val = interpmat * fun(data.inside(:));
  end;
  if (hasmsk && ~strcmp(cfg.projmethod,'project')),
	  % also interpolate the opacity mask
	  maskval = interpmat * msk(data.inside(:));
  end
  
  if ~isempty(cfg.surfinflated)
	  % read the inflated triangulated cortical surface from file
    tmp = load(cfg.surfinflated, 'bnd');
    surf = tmp.bnd;
    if isfield(surf, 'transform'),
      % compute the surface vertices in head coordinates
      surf.pnt = warp_apply(surf.transform, surf.pnt);
    end
  end

  %------do the plotting
  cortex_light = [0.781 0.762 0.664];
  cortex_dark  = [0.781 0.762 0.664]/2;
  if isfield(surf, 'curv')
    % the curvature determines the color of gyri and sulci
    color = surf.curv(:) * cortex_light + (1-surf.curv(:)) * cortex_dark;
  else
    color = repmat(cortex_light, size(surf.pnt,1), 1);
  end

  h1 = patch('Vertices', surf.pnt, 'Faces', surf.tri, 'FaceVertexCData', color , 'FaceColor', 'interp');
  set(h1, 'EdgeColor', 'none');
  axis   off;
  axis vis3d;
  axis equal;

  h2 = patch('Vertices', surf.pnt, 'Faces', surf.tri, 'FaceVertexCData', val , 'FaceColor', 'interp');
  set(h2, 'EdgeColor', 'none');
  if hasmsk
    set(h2, 'FaceVertexAlphaData', maskval);
    set(h2, 'FaceAlpha',          'interp');
    set(h2, 'AlphaDataMapping',   'scaled');
    try
      alim(gca, [opacmin opacmax]);
    end
  end
  try
    caxis(gca,[fcolmin fcolmax]);
  end
  lighting gouraud
  if hasfun
    colormap(cfg.funcolormap);
  end
  if hasmsk
    alphamap(cfg.opacitymap);
  end

  if strcmp(cfg.camlight,'yes')
    camlight
  end

  if strcmp(cfg.colorbar,  'yes'),
    if hasfun
      % use a normal Matlab colorbar
      hc = colorbar;
      set(hc, 'YLim', [fcolmin fcolmax]);
    else
      warning('no colorbar possible without functional data')
    end
  end

elseif isequal(cfg.method,'slice')
  % white BG => mskana

  %% TODO: HERE THE FUNCTION THAT MAKES TO SLICE DIMENSION ALWAYS THE THIRD
  %% DIMENSION, AND ALSO KEEP TRANSFORMATION MATRIX UP TO DATE
  % zoiets
  %if hasana; ana = shiftdim(ana,cfg.slicedim-1); end;
  %if hasfun; fun = shiftdim(fun,cfg.slicedim-1); end;
  %if hasmsk; msk = shiftdim(msk,cfg.slicedim-1); end;
  %%%%% select slices
  if ~ischar(cfg.slicerange)
    ind_fslice = cfg.slicerange(1);
    ind_lslice = cfg.slicerange(2);
  elseif isequal(cfg.slicerange, 'auto')
    if hasfun %default
      if isfield(data,'inside')
        ind_fslice = min(find(max(max(data.inside,[],1),[],2)));
        ind_lslice = max(find(max(max(data.inside,[],1),[],2)));
      else
        ind_fslice = min(find(~isnan(max(max(fun,[],1),[],2))));
        ind_lslice = max(find(~isnan(max(max(fun,[],1),[],2))));
      end
    elseif hasana %if only ana, no fun
      ind_fslice = min(find(max(max(ana,[],1),[],2)));
      ind_lslice = max(find(max(max(ana,[],1),[],2)));
    else
      error('no functional parameter and no anatomical parameter, can not plot');
    end
  else
    error('do not understand cfg.slicerange');
  end
  ind_allslice = linspace(ind_fslice,ind_lslice,cfg.nslices);
  ind_allslice = round(ind_allslice);
  % make new ana, fun, msk, mskana with only the slices that will be plotted (slice dim is always third dimension)
  if hasana; new_ana = ana(:,:,ind_allslice); clear ana; ana=new_ana; clear new_ana; end;
  if hasfun; new_fun = fun(:,:,ind_allslice); clear fun; fun=new_fun; clear new_fun; end;
  if hasmsk; new_msk = msk(:,:,ind_allslice); clear msk; msk=new_msk; clear new_msk; end;
  %if hasmskana; new_mskana = mskana(:,:,ind_allslice); clear mskana; mskana=new_mskana; clear new_mskana; end;

  % update the dimensions of the volume
  if hasana; dim=size(ana); else dim=size(fun); end;

  %%%%% make "quilts", that contain all slices on 2D patched sheet
  % Number of patches along sides of Quilt (M and N)
  % Size (in voxels) of side of patches of Quilt (m and n)
  m = dim(1);
  n = dim(2);
  M = ceil(sqrt(dim(3)));
  N = ceil(sqrt(dim(3)));
  num_patch = N*M;
  if cfg.slicedim~=3
    error('only supported for slicedim=3');
  end
  num_slice = (dim(cfg.slicedim));
  num_empt = num_patch-num_slice;
  % put empty slides on ana, fun, msk, mskana to fill Quilt up
  if hasana; ana(:,:,end+1:num_patch)=0; end;
  if hasfun; fun(:,:,end+1:num_patch)=0; end;
  if hasmsk; msk(:,:,end+1:num_patch)=0; end;
  %if hasmskana; mskana(:,:,end:num_patch)=0; end;
  % put the slices in the quilt
  for iSlice = 1:num_slice
    xbeg = floor((iSlice-1)./M);
    ybeg = mod(iSlice-1, M);
    if hasana
      quilt_ana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(ana(:,:,iSlice));
    end
    if hasfun
      quilt_fun(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(fun(:,:,iSlice));
    end
    if hasmsk
      quilt_msk(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(msk(:,:,iSlice));
    end
    %     if hasmskana
    %       quilt_mskana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(mskana(:,:,iSlice));
    %     end
  end
  % make vols and scales, containes volumes to be plotted (fun, ana, msk) %added ingnie
  if hasana; vols2D{1} = quilt_ana; scales{1} = []; end; % needed when only plotting ana
  if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end;
  if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end;
  plot2D(vols2D, scales, doimage);
  axis off
  if strcmp(cfg.colorbar,  'yes'),
    if hasfun
      % use a normal Matlab coorbar
      hc = colorbar;
      set(hc, 'YLim', [fcolmin fcolmax]);
    else
      warning('no colorbar possible without functional data')
    end
  end

end

title(cfg.title);
set(gcf, 'renderer', cfg.renderer);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION makes an overlay of 3D anatomical, functional and probability
% volumes. The three volumes must be scaled between 0 and 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vols2D] = handle_ortho(vols, indx, slicedir, dim, doimage)

% put 2Dvolumes in fun, ana and msk
if length(vols)>=1 && isempty(vols{1}); hasana=0; else ana=vols{1}; hasana=1; end;
if length(vols)>=2
  if isempty(vols{2}); hasfun=0; else fun=vols{2}; hasfun=1; end;
else hasfun=0; end
if length(vols)>=3
  if isempty(vols{3}); hasmsk=0; else msk=vols{3}; hasmsk=1; end;
else hasmsk=0; end

% select the indices of the intersection
xi = indx(1);
yi = indx(2);
zi = indx(3);
qi = indx(4);
if length(indx)>4,
  qi(2) = indx(5);
else
  qi(2) = 1;
end

% select the slice to plot
if slicedir==1
  yi = 1:dim(2);
  zi = 1:dim(3);
elseif slicedir==2
  xi = 1:dim(1);
  zi = 1:dim(3);
elseif slicedir==3
  xi = 1:dim(1);
  yi = 1:dim(2);
end

% cut out the slice of interest
if hasana; ana = squeeze(ana(xi,yi,zi)); end;
if hasfun; 
  if doimage
    fun = squeeze(fun(xi,yi,zi,:));
  else
    fun = squeeze(fun(xi,yi,zi,qi(1),qi(2)));
  end
end
if hasmsk && length(size(msk))>3
  msk = squeeze(msk(xi,yi,zi,qi(1),qi(2)));
elseif hasmsk  
  msk = squeeze(msk(xi,yi,zi));
end;

%put fun, ana and msk in vols2D
if hasana; vols2D{1} = ana; end;
if hasfun; vols2D{2} = fun; end;
if hasmsk; vols2D{3} = msk; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION plots a two dimensional plot, used in ortho and slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot2D(vols2D, scales, doimage)
cla;
% put 2D volumes in fun, ana and msk
hasana = length(vols2D)>0 && ~isempty(vols2D{1});
hasfun = length(vols2D)>1 && ~isempty(vols2D{2});
hasmsk = length(vols2D)>2 && ~isempty(vols2D{3});

% the transpose is needed for displaying the matrix using the Matlab image() function
if hasana; ana = vols2D{1}'; end;
if hasfun && ~doimage; fun = vols2D{2}'; end;
if hasfun && doimage;  fun = permute(vols2D{2},[2 1 3]); end;
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
  ha = imagesc(ana);
end
hold on

if hasfun

  if doimage
    hf = image(fun);
  else
    hf = imagesc(fun);
    try
      caxis(scales{2});
    end
    % apply the opacity mask to the functional data
    if hasmsk
      % set the opacity
      set(hf, 'AlphaData', msk)
      set(hf, 'AlphaDataMapping', 'scaled')
      try
        alim(scales{3});
      end
    elseif hasana
      set(hf, 'AlphaData', 0.5)
    end

  end
end

axis equal
axis tight
axis xy

