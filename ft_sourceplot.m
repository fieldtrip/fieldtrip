function ft_sourceplot(cfg, functional, anatomical)

% FT_SOURCEPLOT plots functional source reconstruction data on slices or on a
% surface, optionally as an overlay on anatomical MRI data, where statistical data
% can be used to determine the opacity of the mask. Input data comes from
% FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE or statistical values from
% FT_SOURCESTATISTICS.
%
% Use as
%   ft_sourceplot(cfg, data)
% where the input data can contain an anatomical MRI, functional source
% reconstruction results and/or statistical data. If both anatomical and
% functional/statistical data is provided as input, they should be represented or
% interpolated on the same same 3-D grid, e.g. using FT_SOURCEINTERPOLATE.
%
% The slice and ortho visualization plot the data in the input data voxel
% arrangement, i.e. the three ortho views are the 1st, 2nd and 3rd dimension of
% the 3-D data matrix, not of the head coordinate system. The specification of the
% coordinate for slice intersection is specified in head coordinates, i.e.
% relative to the fiducials and in mm or cm. If you want the visualisation to be
% consistent with the head coordinate system, you can reslice the data using
% FT_VOLUMERESLICE.
%
% The configuration should contain:
%   cfg.method        = 'slice',      plots the data on a number of slices in the same plane
%                       'ortho',      plots the data on three orthogonal slices
%                       'surface',    plots the data on a 3D brain surface
%                       'glassbrain', plots a max-projection through the brain
%                       'vertex',     plots the grid points or vertices scaled according to the functional value
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
%   cfg.atlas         = string, filename of atlas to use (default = []) see FT_READ_ATLAS
%                        for ROI masking (see "masking" below) or in "ortho-plotting" mode (see "ortho-plotting" below)
%
% The following parameters can be used for the functional data:
%   cfg.funcolormap   = colormap for functional data, see COLORMAP (default = 'auto')
%                       'auto', depends structure funparameter, or on funcolorlim
%                         - funparameter: only positive values, or funcolorlim:'zeromax' -> 'hot'
%                         - funparameter: only negative values, or funcolorlim:'minzero' -> 'cool'
%                         - funparameter: both pos and neg values, or funcolorlim:'maxabs' -> 'default'
%                         - funcolorlim: [min max] if min & max pos-> 'hot', neg-> 'cool', both-> 'default'
%   cfg.funcolorlim   = color range of the functional data (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(funparameter)) to +max(abs(funparameter))
%                        'zeromax', from 0 to max(funparameter)
%                        'minzero', from min(funparameter) to 0
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
% The following parameters apply for ortho-plotting
%   cfg.location      = location of cut, (default = 'auto')
%                        'auto', 'center' if only anatomy, 'max' if functional data
%                        'min' and 'max' position of min/max funparameter
%                        'center' of the brain
%                        [x y z], coordinates in voxels or head, see cfg.locationcoordinates
%   cfg.locationcoordinates = coordinate system used in cfg.location, 'head' or 'voxel' (default = 'head')
%                              'head', headcoordinates as mm or cm
%                              'voxel', voxelcoordinates as indices
%   cfg.crosshair     = 'yes' or 'no' (default = 'yes')
%   cfg.axis          = 'on' or 'off' (default = 'on')
%   cfg.queryrange    = number, in atlas voxels (default 3)
%
%
% The following parameters apply for slice-plotting
%   cfg.nslices       = number of slices, (default = 20)
%   cfg.slicerange    = range of slices in data, (default = 'auto')
%                       'auto', full range of data
%                       [min max], coordinates of first and last slice in voxels
%   cfg.slicedim      = dimension to slice 1 (x-axis) 2(y-axis) 3(z-axis) (default = 3)
%   cfg.title         = string, title of the figure window
%
% When cfg.method = 'surface', the functional data will be rendered onto a cortical
% mesh (can be an inflated mesh). If the input source data contains a tri-field (i.e.
% a description of a mesh), no interpolation is needed. If the input source data does
% not contain a tri-field, an interpolation is performed onto a specified surface.
% Note that the coordinate system in which the surface is defined should be the same
% as the coordinate system that is represented in source.pos.
%
% The following parameters apply to surface-plotting when an interpolation
% is required
%   cfg.surffile       = string, file that contains the surface (default = 'surface_white_both.mat')
%                        'surface_white_both.mat' contains a triangulation that corresponds with the
%                         SPM anatomical template in MNI coordinates
%   cfg.surfinflated   = string, file that contains the inflated surface (default = [])
%                        may require specifying a point-matching (uninflated) surffile
%   cfg.surfdownsample = number (default = 1, i.e. no downsampling)
%   cfg.projmethod     = projection method, how functional volume data is projected onto surface
%                        'nearest', 'project', 'sphere_avg', 'sphere_weighteddistance'
%   cfg.projvec        = vector (in mm) to allow different projections that
%                        are combined with the method specified in cfg.projcomb
%   cfg.projcomb       = 'mean', 'max', method to combine the different projections
%   cfg.projweight     = vector of weights for the different projections (default = 1)
%   cfg.projthresh     = implements thresholding on the surface level
%                        for example, 0.7 means 70% of maximum
%   cfg.sphereradius   = maximum distance from each voxel to the surface to be
%                        included in the sphere projection methods, expressed in mm
%   cfg.distmat        = precomputed distance matrix (default = [])
%
% The following parameters apply to surface-plotting independent of whether
% an interpolation is required
%   cfg.camlight       = 'yes' or 'no' (default = 'yes')
%   cfg.renderer       = 'painters', 'zbuffer', ' opengl' or 'none' (default = 'opengl')
%                        note that when using opacity the OpenGL renderer is required.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat file on
% disk. This mat files should contain only a single variable corresponding to the
% input structure.
%
% See also FT_SOURCEMOVIE, FT_SOURCEANALYSIS, FT_SOURCEGRANDAVERAGE,
% FT_SOURCESTATISTICS, FT_VOLUMELOOKUP, FT_READ_ATLAS, FT_READ_MRI

% TODO have to be built in:
%   cfg.marker        = [Nx3] array defining N marker positions to display (orig: from sliceinterp)
%   cfg.markersize    = radius of markers (default = 5)
%   cfg.markercolor   = [1x3] marker color in RGB (default = [1 1 1], i.e. white) (orig: from sliceinterp)
%   white background option

% undocumented TODO
%   slice in all directions
%   surface also optimal when inside present
%   come up with a good glass brain projection

% Copyright (C) 2007-2016, Robert Oostenveld, Ingrid Nieuwenhuis
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
ft_preamble loadvar functional anatomical
ft_preamble provenance functional anatomical
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% this is not supported any more as of 26/10/2011
if ischar(functional)
  error('please use cfg.inputfile instead of specifying the input variable as a sting');
end

% ensure that old and unsupported options are not being relied on by the end-user's script
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'funparameter', 'avg.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.pow', 'pow'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.coh', 'coh'});
cfg = ft_checkconfig(cfg, 'renamedval', {'maskparameter', 'avg.mom', 'mom'});
cfg = ft_checkconfig(cfg, 'renamedval', {'location', 'interactive', 'auto'});
% instead of specifying cfg.coordsys, the user should specify the coordsys in the functional data
cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'inputcoordsys', 'coordinates'});
cfg = ft_checkconfig(cfg, 'deprecated', 'coordsys');
% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2837
cfg = ft_checkconfig(cfg, 'renamed', {'viewdim', 'axisratio'});

if isfield(cfg, 'atlas') && ~isempty(cfg.atlas)
  % the atlas lookup requires the specification of the coordsys
  functional     = ft_checkdata(functional, 'datatype', {'volume', 'source'}, 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');
else
  % check if the input functional is valid for this function, a coordsys is not directly needed
  functional     = ft_checkdata(functional, 'datatype', {'volume', 'source'}, 'feedback', 'yes', 'hasunit', 'yes');
end

% set the defaults for all methods
cfg.method        = ft_getopt(cfg, 'method',        'ortho');
cfg.funparameter  = ft_getopt(cfg, 'funparameter',  []);
cfg.maskparameter = ft_getopt(cfg, 'maskparameter', []);
cfg.downsample    = ft_getopt(cfg, 'downsample',    1);
cfg.title         = ft_getopt(cfg, 'title',         '');
cfg.atlas         = ft_getopt(cfg, 'atlas',         []);
cfg.marker        = ft_getopt(cfg, 'marker',        []);
cfg.markersize    = ft_getopt(cfg, 'markersize',    5);
cfg.markercolor   = ft_getopt(cfg, 'markercolor',   [1 1 1]);
cfg.renderer      = ft_getopt(cfg, 'renderer',      'opengl');
cfg.colorbar      = ft_getopt(cfg, 'colorbar',      'yes');
cfg.voxelratio    = ft_getopt(cfg, 'voxelratio',    'data'); % display size of the voxel, 'data' or 'square'
cfg.axisratio     = ft_getopt(cfg, 'axisratio',     'data'); % size of the axes of the three orthoplots, 'square', 'voxel', or 'data'

if ~isfield(cfg, 'anaparameter')
  if isfield(functional, 'anatomy')
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

if isfield(cfg, 'TTlookup'),
  error('TTlookup is old; now specify cfg.atlas, see help!');
end

% select the functional and the mask parameter
cfg.funparameter  = parameterselection(cfg.funparameter, functional);
cfg.maskparameter = parameterselection(cfg.maskparameter, functional);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1};  end
try, cfg.maskparameter = cfg.maskparameter{1}; end


% the data can be passed as input argument or can be read from disk
hasanatomical = exist('anatomical', 'var');

if hasanatomical
  % interpolate on the fly, this also does the downsampling if requested
  tmpcfg = keepfields(cfg, {'downsample', 'interpmethod'});
  tmpcfg.parameter = cfg.funparameter;
  functional = ft_sourceinterpolate(tmpcfg, functional, anatomical);
  [cfg, functional] = rollback_provenance(cfg, functional);
elseif ~hasanatomical && cfg.downsample~=1
  % optionally downsample the functional volume
  tmpcfg = keepfields(cfg, {'downsample'});
  tmpcfg.parameter = {cfg.funparameter, cfg.maskparameter, cfg.anaparameter};
  functional = ft_volumedownsample(tmpcfg, functional);
  [cfg, functional] = rollback_provenance(cfg, functional);
end

if isfield(functional, 'dim') && isfield(functional, 'transform')
  % this is a regular 3D functional volume
  isUnstructuredFun = false;
elseif isfield(functional, 'dim') && isfield(functional, 'pos')
  % these are positions that can be mapped onto a 3D regular grid
  isUnstructuredFun  = false;
  % contstruct the transformation matrix from the positions
  functional.transform = pos2transform(functional.pos, functional.dim);
else
  % this is functional data on irregular positions, such as a cortical sheet
  isUnstructuredFun = true;
end

% this only relates to the dimensions of the geometry, which is npos*1 or nx*ny*nz
if isUnstructuredFun
  dim = [size(functional.pos,1) 1];
else
  dim = functional.dim;
end


%% get the elements that will be plotted

hasatlas = ~isempty(cfg.atlas);
if hasatlas
  if ischar(cfg.atlas)
    % initialize the atlas
    [p, f, x] = fileparts(cfg.atlas);
    fprintf(['reading ', f, ' atlas coordinates and labels\n']);
    atlas = ft_read_atlas(cfg.atlas);
  else
    atlas = cfg.atlas;
  end
end

hasroi = ~isempty(cfg.roi);
if hasroi
  if ~hasatlas
    error('specify cfg.atlas which specifies cfg.roi')
  else
    % get the mask
    tmpcfg          = [];
    tmpcfg.roi      = cfg.roi;
    tmpcfg.atlas    = cfg.atlas;
    tmpcfg.inputcoord = functional.coordsys;
    roi = ft_volumelookup(tmpcfg,functional);
  end
end

hasana = isfield(functional, cfg.anaparameter);
if hasana
  ana = getsubfield(functional, cfg.anaparameter);
  if isa(ana, 'uint8') || isa(ana, 'uint16') || isa(ana, 'int8') || isa(ana, 'int16')
    ana = double(ana);
  end
  fprintf('scaling anatomy to [0 1]\n');
  dmin = min(ana(:));
  dmax = max(ana(:));
  ana  = (ana-dmin)./(dmax-dmin);
  ana  = reshape(ana, dim);
end

%%% funparameter
hasfun = isfield(functional, cfg.funparameter);
if hasfun
  fun = getsubfield(functional, cfg.funparameter);

  dimord = getdimord(functional, cfg.funparameter);
  dimtok = tokenize(dimord, '_');

  % replace the cell-array functional with a normal array
  if strcmp(dimtok{1}, '{pos}')
    tmpdim = getdimsiz(functional, cfg.funparameter);
    tmpfun = nan(tmpdim);
    insideindx = find(functional.inside);
    for i=insideindx(:)'
      tmpfun(i,:) = fun{i};
    end
    fun = tmpfun;
    clear tmpfun
    dimtok{1} = 'pos';  % update the description of the dimensions
    dimord([1 5]) = []; % remove the { and }
  end

  % ensure that the functional data is real
  if ~isreal(fun)
    warning('functional data is complex, taking absolute value');
    fun = abs(fun);
  end

  if strcmp(dimord, 'pos_rgb')
    % treat functional data as rgb values
    if any(fun(:)>1 | fun(:)<0)
      % scale
      tmpdim = size(fun);
      nvox   = prod(tmpdim(1:end-1));
      tmpfun = reshape(fun,[nvox tmpdim(end)]);
      m1     = max(tmpfun,[],1);
      m2     = min(tmpfun,[],1);
      tmpfun = (tmpfun-m2(ones(nvox,1),:))./(m1(ones(nvox,1),:)-m2(ones(nvox,1),:));
      fun    = reshape(tmpfun, tmpdim);
      clear tmpfun
    end
    qi      = 1;
    hasfreq = 0;
    hastime = 0;

    doimage = 1;
    fcolmin = 0;
    fcolmax = 1;

  else
    % determine scaling min and max (fcolmin fcolmax) and funcolormap
    if ~isa(fun, 'logical')
      funmin = min(fun(:));
      funmax = max(fun(:));
    else
      funmin = 0;
      funmax = 1;
    end
    % smart automatic limits
    if isequal(cfg.funcolorlim, 'auto')
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
      if isequal(cfg.funcolorlim, 'maxabs')
        fcolmin = -max(abs([funmin,funmax]));
        fcolmax =  max(abs([funmin,funmax]));
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'default'; end;
      elseif isequal(cfg.funcolorlim, 'zeromax')
        fcolmin = 0;
        fcolmax = funmax;
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'hot'; end;
      elseif isequal(cfg.funcolorlim, 'minzero')
        fcolmin = funmin;
        fcolmax = 0;
        if isequal(cfg.funcolormap, 'auto'); cfg.funcolormap = 'cool'; end;
      else
        error('do not understand cfg.funcolorlim');
      end
    else
      % limits are numeric
      fcolmin = cfg.funcolorlim(1);
      fcolmax = cfg.funcolorlim(2);
      % smart colormap
      if isequal(cfg.funcolormap, 'auto')
        if sign(fcolmin) == -1 && sign(fcolmax) == 1
          cfg.funcolormap = 'default';
        else
          if fcolmin < 0
            cfg.funcolormap = 'cool';
          else
            cfg.funcolormap = 'hot';
          end
        end
      end
    end % if ischar
    clear funmin funmax

    % what if fun is 4D?
    if ndims(fun)>3 || prod(dim)==size(fun,1)
      if strcmp(dimord, 'pos_freq_time')
        % functional contains time-frequency representation
        qi      = [1 1];
        hasfreq = numel(functional.freq)>1;
        hastime = numel(functional.time)>1;
        fun     = reshape(fun, [dim numel(functional.freq) numel(functional.time)]);
      elseif strcmp(dimord, 'pos_time')
        % functional contains evoked field
        qi      = 1;
        hasfreq = 0;
        hastime = numel(functional.time)>1;
        fun     = reshape(fun, [dim numel(functional.time)]);
      elseif strcmp(dimord, 'pos_freq')
        % functional contains frequency spectra
        qi      = 1;
        hasfreq = numel(functional.freq)>1;
        hastime = 0;
        fun     = reshape(fun, [dim numel(functional.freq)]);
      else
        qi      = 1;
        hasfreq = 0;
        hastime = 0;
        fun     = reshape(fun, dim);
      end
    else
      % do nothing
      qi      = 1;
      hasfreq = 0;
      hastime = 0;
    end

    doimage = 0;
  end % if dimord has rgb or something else

else
  % there is no functional data
  qi      = 1;
  hasfreq = 0;
  hastime = 0;
  doimage = 0;
  fcolmin = 0; % needs to be defined for callback
  fcolmax = 1;
end

hasmsk = issubfield(functional, cfg.maskparameter);
if hasmsk
  if ~hasfun
    error('you can not have a mask without functional data')
  else
    msk = getsubfield(functional, cfg.maskparameter);
    if islogical(msk) % otherwise sign() not posible
      msk = double(msk);
    end
  end
  % reshape to match fun
  if strcmp(dimord, 'pos_freq_time')
    % functional contains timefrequency representation
    msk     = reshape(msk, [dim numel(functional.freq) numel(functional.time)]);
  elseif strcmp(dimord, 'pos_time')
    % functional contains evoked field
    msk     = reshape(msk, [dim numel(functional.time)]);
  elseif strcmp(dimord, 'pos_freq')
    % functional contains frequency spectra
    msk     = reshape(msk, [dim numel(functional.freq)]);
  else
    msk     = reshape(msk, dim);
  end

  % determine scaling and opacitymap
  mskmin = min(msk(:));
  mskmax = max(msk(:));
  % determine the opacity limits and the opacity map
  % smart limits: make from auto other string, or equal to funcolorlim if funparameter == maskparameter
  if isequal(cfg.opacitylim, 'auto')
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
        if isequal(cfg.opacitymap, 'auto'), cfg.opacitymap = 'rampup'; end;
      case 'minzero'
        opacmin = mskmin;
        opacmax = 0;
        if isequal(cfg.opacitymap, 'auto'), cfg.opacitymap = 'rampdown'; end;
      case 'maxabs'
        opacmin = -max(abs([mskmin, mskmax]));
        opacmax =  max(abs([mskmin, mskmax]));
        if isequal(cfg.opacitymap, 'auto'), cfg.opacitymap = 'vdown'; end;
      otherwise
        error('incorrect specification of cfg.opacitylim');
    end % switch opacitylim
  else
    % limits are numeric
    opacmin = cfg.opacitylim(1);
    opacmax = cfg.opacitylim(2);
    if isequal(cfg.opacitymap, 'auto')
      if sign(opacmin)>-1 && sign(opacmax)>-1
        cfg.opacitymap = 'rampup';
      elseif sign(opacmin)<1 && sign(opacmax)<1
        cfg.opacitymap = 'rampdown';
      else
        cfg.opacitymap = 'vdown';
      end
    end
  end % handling opacitylim and opacitymap
  clear mskmin mskmax
else
  opacmin = [];
  opacmax = [];
end

% prevent outside fun from being plotted
if hasfun && ~hasmsk && isfield(functional, 'inside')
  hasmsk = 1;
  msk = zeros(dim);
  cfg.opacitymap = 'rampup';
  opacmin = 0;
  opacmax = 1;
  % make intelligent mask
  if isequal(cfg.method, 'surface')
    msk(functional.inside) = 1;
  else
    if hasana
      msk(functional.inside) = 0.5; % so anatomy is visible
    else
      msk(functional.inside) = 1;
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
  opacmin = [];
  opacmax = []; % has to be defined
elseif hasroi
  error('you can not have a roi without functional data')
end

%% give some feedback

if ~hasfun && ~hasana
  % this seems to be a problem that people often have due to incorrect specification of the cfg
  error('no anatomy is present and no functional data is selected, please check your cfg.funparameter');
end
if ~hasana
  fprintf('not plotting anatomy\n');
end
if ~hasfun
  fprintf('not plotting functional data\n');
end
if ~hasmsk
  fprintf('not applying a mask on the functional data\n');
end
if ~hasatlas
  fprintf('not using an atlas\n');
end
if ~hasroi
  fprintf('not using a region-of-interest\n');
end


%% start building the figure
h = figure;
set(h, 'color', [1 1 1]);
set(h, 'visible', 'on');
set(h, 'renderer', cfg.renderer);
if ~isempty(cfg.title)
  title(cfg.title);
end

%%% set color and opacity mapping for this figure
if hasfun
  colormap(cfg.funcolormap);
  cfg.funcolormap = colormap;
end
if hasmsk
  cfg.opacitymap = alphamap(cfg.opacitymap);
  alphamap(cfg.opacitymap);
  if ndims(fun)>3 && ndims(msk)==3
    siz = size(fun);
    msk = repmat(msk, [1 1 1 siz(4:end)]);
  end
end

switch cfg.method
  case 'slice'
    % set the defaults for method=slice
    cfg.nslices    = ft_getopt(cfg, 'nslices',    20);
    cfg.slicedim   = ft_getopt(cfg, 'slicedim',   3);
    cfg.slicerange = ft_getopt(cfg, 'slicerange', 'auto');



    % white BG => mskana

    % TODO: HERE THE FUNCTION THAT MAKES TO SLICE DIMENSION ALWAYS THE THIRD DIMENSION, AND ALSO KEEP TRANSFORMATION MATRIX UP TO DATE
    % zoiets
    % if hasana; ana = shiftdim(ana,cfg.slicedim-1); end;
    % if hasfun; fun = shiftdim(fun,cfg.slicedim-1); end;
    % if hasmsk; msk = shiftdim(msk,cfg.slicedim-1); end;

    % ADDED BY JM: allow for slicedim different than 3
    switch cfg.slicedim
      case 1
        dim = dim([2 3 1]);
        if hasana, ana = permute(ana,[2 3 1]); end
        if hasfun, fun = permute(fun,[2 3 1]); end
        if hasmsk, msk = permute(msk,[2 3 1]); end
        cfg.slicedim=3;
      case 2
        dim = dim([3 1 2]);
        if hasana, ana = permute(ana,[3 1 2]); end
        if hasfun, fun = permute(fun,[3 1 2]); end
        if hasmsk, msk = permute(msk,[3 1 2]); end
        cfg.slicedim=3;
      otherwise
        % nothing needed
    end

    %%%%% select slices
    if ~ischar(cfg.slicerange)
      ind_fslice = cfg.slicerange(1);
      ind_lslice = cfg.slicerange(2);
    elseif isequal(cfg.slicerange, 'auto')
      if hasfun % default
        if isfield(functional, 'inside')

          insideMask = false(size(fun));
          insideMask(functional.inside) = true;

          ind_fslice = min(find(max(max(insideMask,[],1),[],2)));
          ind_lslice = max(find(max(max(insideMask,[],1),[],2)));
        else
          ind_fslice = min(find(~isnan(max(max(fun,[],1),[],2))));
          ind_lslice = max(find(~isnan(max(max(fun,[],1),[],2))));
        end
      elseif hasana % if only ana, no fun
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
    % if hasmskana; new_mskana = mskana(:,:,ind_allslice); clear mskana; mskana=new_mskana; clear new_mskana; end;

    % update the dimensions of the volume
    if hasana; dim=size(ana); else dim=size(fun); end;

    %%%%% make a "quilt", that contain all slices on 2D patched sheet
    % Number of patches along sides of Quilt (M and N)
    % Size (in voxels) of side of patches of Quilt (m and n)

    % take care of a potential singleton 3rd dimension
    if numel(dim)<3
      dim(end+1:3) = 1;
    end

    %if cfg.slicedim~=3
    %  error('only supported for slicedim=3');
    %end


    m = dim(1);
    n = dim(2);
    M = ceil(sqrt(dim(3)));
    N = ceil(sqrt(dim(3)));
    num_patch = N*M;

    num_slice = (dim(cfg.slicedim));
    num_empt = num_patch-num_slice;
    % put empty slides on ana, fun, msk, mskana to fill Quilt up
    if hasana; ana(:,:,end+1:num_patch)=0; end;
    if hasfun; fun(:,:,end+1:num_patch)=0; end;
    if hasmsk; msk(:,:,end+1:num_patch)=0; end;
    % if hasmskana; mskana(:,:,end:num_patch)=0; end;
    % put the slices in the quilt
    for iSlice = 1:num_slice
      xbeg = floor((iSlice-1)./M);
      ybeg = mod(iSlice-1, M);
      if hasana
        quilt_ana(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=squeeze(ana(:,:,iSlice));
      end
      if hasfun
        quilt_fun(ybeg*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=squeeze(fun(:,:,iSlice));
      end
      if hasmsk
        quilt_msk(ybeg.*m+1:(ybeg+1)*m, xbeg*n+1:(xbeg+1)*n)=squeeze(msk(:,:,iSlice));
      end
      %     if hasmskana
      %       quilt_mskana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(mskana(:,:,iSlice));
      %     end
    end
    % make vols and scales, containes volumes to be plotted (fun, ana, msk), added by ingnie
    if hasana; vols2D{1} = quilt_ana; scales{1} = []; end; % needed when only plotting ana
    if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end;
    if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end;

    % the transpose is needed for displaying the matrix using the MATLAB image() function
    if hasana;             ana = vols2D{1}'; end;
    if hasfun && ~doimage; fun = vols2D{2}'; end;
    if hasfun &&  doimage; fun = permute(vols2D{2},[2 1 3]); end;
    if hasmsk;             msk = vols2D{3}'; end;

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
    axis off

    if istrue(cfg.colorbar)
      if hasfun
        % use a normal MATLAB coorbar
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
      else
        warning('no colorbar possible without functional data')
      end
    end

  case 'ortho'
    % set the defaults for method=ortho
    cfg.location            = ft_getopt(cfg, 'location',            'auto');
    cfg.locationcoordinates = ft_getopt(cfg, 'locationcoordinates', 'head');
    cfg.crosshair           = ft_getopt(cfg, 'crosshair',           'yes');
    cfg.axis                = ft_getopt(cfg, 'axis',                'on');
    cfg.queryrange          = ft_getopt(cfg, 'queryrange',          3);

    if ~ischar(cfg.location)
      if strcmp(cfg.locationcoordinates, 'head')
        % convert the headcoordinates location into voxel coordinates
        loc = inv(functional.transform) * [cfg.location(:); 1];
        loc = round(loc(1:3));
      elseif strcmp(cfg.locationcoordinates, 'voxel')
        % the location is already in voxel coordinates
        loc = round(cfg.location(1:3));
      else
        error('you should specify cfg.locationcoordinates');
      end
    else
      if isequal(cfg.location, 'auto')
        if hasfun
          if isequal(cfg.funcolorlim, 'maxabs');
            loc = 'max';
          elseif isequal(cfg.funcolorlim, 'zeromax');
            loc = 'max';
          elseif isequal(cfg.funcolorlim, 'minzero');
            loc = 'min';
          else % if numerical
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
      [dummy, minindx] = min(fun(:));
      [xi, yi, zi] = ind2sub(dim, minindx);
    elseif ischar(loc) && strcmp(loc, 'max')
      if isempty(cfg.funparameter)
        error('cfg.location is max, but no functional parameter specified');
      end
      [dummy, maxindx] = max(fun(:));
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

    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

    % axes settings
    if strcmp(cfg.axisratio, 'voxel')
      % determine the number of voxels to be plotted along each axis
      axlen1 = dim(1);
      axlen2 = dim(2);
      axlen3 = dim(3);
    elseif strcmp(cfg.axisratio, 'data')
      % determine the length of the edges along each axis
      [cp_voxel, cp_head] = cornerpoints(dim, functional.transform);
      axlen1 = norm(cp_head(2,:)-cp_head(1,:));
      axlen2 = norm(cp_head(4,:)-cp_head(1,:));
      axlen3 = norm(cp_head(5,:)-cp_head(1,:));
    elseif strcmp(cfg.axisratio, 'square')
      % the length of the axes should be equal
      axlen1 = 1;
      axlen2 = 1;
      axlen3 = 1;
    end

    % this is the size reserved for subplot h1, h2 and h3
    h1size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h1size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h2size(1) = 0.82*axlen2/(axlen1 + axlen2);
    h2size(2) = 0.82*axlen3/(axlen2 + axlen3);
    h3size(1) = 0.82*axlen1/(axlen1 + axlen2);
    h3size(2) = 0.82*axlen2/(axlen2 + axlen3);

    if strcmp(cfg.voxelratio, 'square')
      voxlen1 = 1;
      voxlen2 = 1;
      voxlen3 = 1;
    elseif strcmp(cfg.voxelratio, 'data')
      % the size of the voxel is scaled with the data
      [cp_voxel, cp_head] = cornerpoints(dim, functional.transform);
      voxlen1 = norm(cp_head(2,:)-cp_head(1,:))/norm(cp_voxel(2,:)-cp_voxel(1,:));
      voxlen2 = norm(cp_head(4,:)-cp_head(1,:))/norm(cp_voxel(4,:)-cp_voxel(1,:));
      voxlen3 = norm(cp_head(5,:)-cp_head(1,:))/norm(cp_voxel(5,:)-cp_voxel(1,:));
    end

    %% the figure is interactive, add callbacks
    set(h, 'windowbuttondownfcn', @cb_buttonpress);
    set(h, 'windowbuttonupfcn',   @cb_buttonrelease);
    set(h, 'windowkeypressfcn',   @cb_keyboard);
    set(h, 'CloseRequestFcn',     @cb_cleanup);

    % ensure that this is done in interactive mode
    set(h, 'renderer', cfg.renderer);

    %% create figure handles

    % axis handles will hold the anatomical functional if present, along with labels etc.
    h1 = axes('position',[0.06                0.06+0.06+h3size(2) h1size(1) h1size(2)]);
    h2 = axes('position',[0.06+0.06+h1size(1) 0.06+0.06+h3size(2) h2size(1) h2size(2)]);
    h3 = axes('position',[0.06                0.06                h3size(1) h3size(2)]);

    set(h1, 'Tag', 'ik', 'Visible', cfg.axis, 'XAxisLocation', 'top');
    set(h2, 'Tag', 'jk', 'Visible', cfg.axis, 'YAxisLocation', 'right'); % after rotating in ft_plot_ortho this becomes top
    set(h3, 'Tag', 'ij', 'Visible', cfg.axis);

    set(h1, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    set(h2, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);
    set(h3, 'DataAspectRatio',1./[voxlen1 voxlen2 voxlen3]);

    % create structure to be passed to gui
    opt               = [];
    opt.dim           = dim;
    opt.ijk           = [xi yi zi];
    opt.h1size        = h1size;
    opt.h2size        = h2size;
    opt.h3size        = h3size;
    opt.handlesaxes   = [h1 h2 h3];
    opt.handlesfigure = h;
    opt.axis          = cfg.axis;
    if hasatlas
      opt.atlas = atlas;
    end
    if hasana
      opt.ana = ana;
    end
    if hasfun
      opt.fun = fun;
    end
    if hasmsk
      opt.msk = msk;
    end
    opt.update        = [1 1 1];
    opt.init          = true;
    opt.usedim        = (isUnstructuredFun==false);
    opt.usepos        = (isUnstructuredFun==true);
    opt.hasatlas      = hasatlas;
    opt.hasfreq       = hasfreq;
    opt.hastime       = hastime;
    opt.hasmsk        = hasmsk;
    opt.hasfun        = hasfun;
    opt.hasana        = hasana;
    opt.qi            = qi;
    opt.tag           = 'ik';
    opt.functional    = functional;
    opt.fcolmin       = fcolmin;
    opt.fcolmax       = fcolmax;
    opt.opacmin       = opacmin;
    opt.opacmax       = opacmax;
    opt.clim          = []; % contrast limits for the anatomy, see ft_volumenormalise
    opt.colorbar      = cfg.colorbar;
    opt.queryrange    = cfg.queryrange;
    opt.funcolormap   = cfg.funcolormap;
    opt.crosshair     = istrue(cfg.crosshair);

    %% do the actual plotting
    setappdata(h, 'opt', opt);
    cb_redraw(h);

    fprintf('\n');
    fprintf('click left mouse button to reposition the cursor\n');
    fprintf('click and hold right mouse button to update the position while moving the mouse\n');
    fprintf('use the arrowkeys to navigate in the current axis\n');


  case 'surface'
    % set the defaults for method=surface
    cfg.downsample     = ft_getopt(cfg, 'downsample',     1);
    cfg.surfdownsample = ft_getopt(cfg, 'surfdownsample', 1);
    cfg.surffile       = ft_getopt(cfg, 'surffile', 'surface_white_both.mat'); % use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
    cfg.surfinflated   = ft_getopt(cfg, 'surfinflated',  []);
    cfg.sphereradius   = ft_getopt(cfg, 'sphereradius',  []);
    cfg.projvec        = ft_getopt(cfg, 'projvec',       1);
    cfg.projweight     = ft_getopt(cfg, 'projweight',    ones(size(cfg.projvec)));
    cfg.projcomb       = ft_getopt(cfg, 'projcomb',      'mean'); % or max
    cfg.projthresh     = ft_getopt(cfg, 'projthresh',    []);
    cfg.projmethod     = ft_getopt(cfg, 'projmethod',    'nearest');
    cfg.distmat        = ft_getopt(cfg, 'distmat',       []);
    cfg.camlight       = ft_getopt(cfg, 'camlight',      'yes');

    % determine whether the source functional already contains a triangulation
    interpolate2surf = 0;
    if ~isUnstructuredFun
      % no triangulation present: interpolation should be performed
      fprintf('The source functional is defined on a 3D grid, interpolation to a surface mesh will be performed\n');
      interpolate2surf = 1;
    elseif isUnstructuredFun && isfield(functional, 'tri')
      fprintf('The source functional is defined on a triangulated surface, using the surface mesh description in the functional\n');
    elseif isUnstructuredFun
      % add a transform field to the functional
      fprintf('The source functional does not contain a triangulated surface, we may need to interpolate to a surface mesh\n');
      functional.transform = pos2transform(functional.pos);
      interpolate2surf = 1;
    end

    if interpolate2surf,
      % deal with the interpolation
      % FIXME this should be dealt with by ft_sourceinterpolate

      % read the triangulated cortical surface from file
      surf = ft_read_headshape(cfg.surffile);

      if isfield(surf, 'transform'),
        % compute the surface vertices in head coordinates
        surf.pos = ft_warp_apply(surf.transform, surf.pos);
      end

      % downsample the cortical surface
      if cfg.surfdownsample > 1
        if ~isempty(cfg.surfinflated)
          error('downsampling the surface is not possible in combination with an inflated surface');
        end
        fprintf('downsampling surface from %d vertices\n', size(surf.pos,1));
        [temp.tri, temp.pos] = reducepatch(surf.tri, surf.pos, 1/cfg.surfdownsample);
        % find indices of retained patch faces
        [dummy, idx] = ismember(temp.pos, surf.pos, 'rows');
        idx(idx==0)  = [];
        surf.tri = temp.tri;
        surf.pos = temp.pos;
        clear temp
        % downsample other fields
        if isfield(surf, 'curv'),       surf.curv       = surf.curv(idx);       end
        if isfield(surf, 'sulc'),       surf.sulc       = surf.sulc(idx);       end
        if isfield(surf, 'hemisphere'), surf.hemisphere = surf.hemisphere(idx); end
      end

      % these are required
      if ~isfield(functional, 'inside')
        functional.inside = true(dim);
      end

      fprintf('%d voxels in functional data\n', prod(dim));
      fprintf('%d vertices in cortical surface\n', size(surf.pos,1));

      tmpcfg = [];
      tmpcfg.parameter = {cfg.funparameter};
      if ~isempty(cfg.maskparameter)
        tmpcfg.parameter = [tmpcfg.parameter {cfg.maskparameter}];
        maskparameter    = cfg.maskparameter;
      else
        tmpcfg.parameter = [tmpcfg.parameter {'mask'}];
        functional.mask  = msk;
        maskparameter    = 'mask'; % temporary variable
      end
      tmpcfg.interpmethod = cfg.projmethod;
      tmpcfg.distmat      = cfg.distmat;
      tmpcfg.sphereradius = cfg.sphereradius;
      tmpcfg.projvec      = cfg.projvec;
      tmpcfg.projcomb     = cfg.projcomb;
      tmpcfg.projweight   = cfg.projweight;
      tmpcfg.projthresh   = cfg.projthresh;
      tmpdata             = ft_sourceinterpolate(tmpcfg, functional, surf);

      if hasfun, val      = getsubfield(tmpdata, cfg.funparameter);  val     = val(:);     end
      if hasmsk, maskval  = getsubfield(tmpdata, maskparameter);     maskval = maskval(:); end

      if ~isempty(cfg.projthresh),
        maskval(abs(val) < cfg.projthresh*max(abs(val(:)))) = 0;
      end

    else
      surf     = [];
      surf.pos = functional.pos;
      surf.tri = functional.tri;

      % if hasfun, val     = fun(functional.inside(:)); end
      % if hasmsk, maskval = msk(functional.inside(:)); end
      if hasfun, val     = fun(:); end
      if hasmsk, maskval = msk(:); end

    end

    if ~isempty(cfg.surfinflated)
      if ~isstruct(cfg.surfinflated)
        % read the inflated triangulated cortical surface from file
        surf = ft_read_headshape(cfg.surfinflated);
      else
        surf = cfg.surfinflated;
        if isfield(surf, 'transform'),
          % compute the surface vertices in head coordinates
          surf.pos = ft_warp_apply(surf.transform, surf.pos);
        end
      end
    end

    %------do the plotting
    cortex_light = [0.781 0.762 0.664];
    cortex_dark  = [0.781 0.762 0.664]/2;
    if isfield(surf, 'curv')
      % the curvature determines the color of gyri and sulci
      color = surf.curv(:) * cortex_dark + (1-surf.curv(:)) * cortex_light;
    else
      color = repmat(cortex_light, size(surf.pos,1), 1);
    end

    h1 = patch('Vertices', surf.pos, 'Faces', surf.tri, 'FaceVertexCData', color, 'FaceColor', 'interp');
    set(h1, 'EdgeColor', 'none');
    axis   off;
    axis vis3d;
    axis equal;

    if hasfun
      h2 = patch('Vertices', surf.pos, 'Faces', surf.tri, 'FaceVertexCData', val, 'FaceColor', 'interp');
      set(h2, 'EdgeColor', 'none');
      try
        caxis(gca,[fcolmin fcolmax]);
      end
      colormap(cfg.funcolormap);
      if hasmsk
        set(h2, 'FaceVertexAlphaData', maskval);
        set(h2, 'FaceAlpha',          'interp');
        set(h2, 'AlphaDataMapping',   'scaled');
        try
          alim(gca, [opacmin opacmax]);
        end
        alphamap(cfg.opacitymap);
      end
    end

    lighting gouraud

    if istrue(cfg.camlight)
      camlight
    end

    if istrue(cfg.colorbar)
      if hasfun
        % use a normal MATLAB colorbar
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
      else
        warning('no colorbar possible without functional data')
      end
    end

  case 'glassbrain'
    % This is implemented using a recursive call with an updated functional data
    % structure. The functional volume is replaced by a volume in which the maxima
    % are projected to the "edge" of the volume.
    tmpcfg                      = keepfields(cfg, {'funparameter', 'funcolorlim', 'funcolormap', 'opacitylim', 'axis', 'renderer'});
    tmpcfg.method               = 'ortho';
    tmpcfg.location             = [1 1 1];
    tmpcfg.locationcoordinates  = 'voxel';
    tmpcfg.maskparameter        = 'inside';

    if hasfun,
      fun = getsubfield(functional, cfg.funparameter);
      fun = reshape(fun, dim);
      fun(1,:,:) = max(fun, [], 1); % get the projection along the 1st dimension
      fun(:,1,:) = max(fun, [], 2); % get the projection along the 2nd dimension
      fun(:,:,1) = max(fun, [], 3); % get the projection along the 3rd dimension
      functional = setsubfield(functional, cfg.funparameter, fun);
    end

    if hasana,
      ana = getsubfield(functional, cfg.anaparameter);
      % this remains as it is
      functional = setsubfield(functional, cfg.anaparameter, ana);
    end

    if hasmsk,
      msk = getsubfield(functional, 'inside');
      msk = reshape(msk, dim);
      if hasana
        msk(1,:,:) = squeeze(fun(1,:,:))>0 & imfill(abs(squeeze(ana(1,:,:))-1))>0;
        msk(:,1,:) = squeeze(fun(:,1,:))>0 & imfill(abs(squeeze(ana(:,1,:))-1))>0;
        msk(:,:,1) = squeeze(fun(:,:,1))>0 & imfill(abs(ana(:,:,1)-1))>0;
      else
        msk(1,:,:) = squeeze(fun(1,:,:))>0;
        msk(:,1,:) = squeeze(fun(:,1,:))>0;
        msk(:,:,1) = squeeze(fun(:,:,1))>0;
      end
      functional = setsubfield(functional, 'inside', msk);
    end

    ft_sourceplot(tmpcfg, functional);

  case 'vertex'
    if isUnstructuredFun
      pos = functional.pos;
    else
      [X, Y, Z] = ndgrid(1:dim(1), 1:dim(2), 1:dim(3));
      pos = ft_warp_apply(functional.transform, [X(:) Y(:) Z(:)]);
    end

    if isfield(functional, 'inside')
      pos = pos(functional.inside,:);
      if hasfun
        fun = fun(functional.inside);
      end
    end

    % scale the functional data between -30 and 30
    fun = 30*fun/max(abs(fun(:)));
    if any(fun<=0)
      warning('using red for positive and blue for negative functional values')
      col = zeros(numel(fun), 3); % RGB
      col(fun>0,1) = 1;  % red
      col(fun<0,3) = 1;  % blue
      fun(fun==0) = eps; % these will be black
      ft_plot_mesh(pos, 'vertexsize', abs(fun), 'vertexcolor', col);
    else
      ft_plot_mesh(pos, 'vertexsize', fun, 'vertexcolor', 'k');
    end

    % ensure that the axes don't change if you rotate
    axis vis3d

  otherwise
    error('unsupported method "%s"', cfg.method);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous functional
ft_postamble provenance

% add a menu to the figure
% also, delete any possibly existing previous menu, this is safe because delete([]) does nothing
ftmenu = uimenu(gcf, 'Label', 'FieldTrip');
uimenu(ftmenu, 'Label', 'Show pipeline',  'Callback', {@menu_pipeline, cfg});
uimenu(ftmenu, 'Label', 'About',  'Callback', @menu_about);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_redraw(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
tag = get(curr_ax, 'tag');

functional = opt.functional;

h1 = opt.handlesaxes(1);
h2 = opt.handlesaxes(2);
h3 = opt.handlesaxes(3);

xi = opt.ijk(1);
yi = opt.ijk(2);
zi = opt.ijk(3);
qi = opt.qi;

if any([xi yi zi] > functional.dim) || any([xi yi zi] <= 0)
  return;
end

opt.ijk = [xi yi zi 1]';
if opt.usedim
  xyz = functional.transform * opt.ijk;
elseif opt.usepos
  ix  = sub2ind(opt.dim,xi,yi,zi);
  xyz = functional.pos(ix,:);
end
opt.ijk = opt.ijk(1:3);

% construct a string with user feedback
str1 = sprintf('voxel %d, indices [%d %d %d]', sub2ind(functional.dim(1:3), xi, yi, zi), opt.ijk);

if isfield(functional, 'coordsys') && isfield(functional, 'unit')
  str2 = sprintf('%s coordinates [%.1f %.1f %.1f] %s', functional.coordsys, xyz(1:3), functional.unit);
elseif ~isfield(functional, 'coordsys') && isfield(functional, 'unit')
  str2 = sprintf('location [%.1f %.1f %.1f] %s', xyz(1:3), functional.unit);
elseif isfield(functional, 'coordsys') && ~isfield(functional, 'unit')
  str2 = sprintf('%s coordinates [%.1f %.1f %.1f]', functional.coordsys, xyz(1:3));
elseif ~isfield(functional, 'coordsys') && ~isfield(functional, 'unit')
  str2 = sprintf('location [%.1f %.1f %.1f]', xyz(1:3));
else
  str2 = '';
end

if opt.hasfreq && opt.hastime,
  str3 = sprintf('%.1f s, %.1f Hz', functional.time(opt.qi(2)), functional.freq(opt.qi(1)));
elseif ~opt.hasfreq && opt.hastime,
  str3 = sprintf('%.1f s', functional.time(opt.qi(1)));
elseif opt.hasfreq && ~opt.hastime,
  str3 = sprintf('%.1f Hz', functional.freq(opt.qi(1)));
else
  str3 = '';
end

if opt.hasfun
  if ~opt.hasfreq && ~opt.hastime
    val = opt.fun(xi, yi, zi);
  elseif ~opt.hasfreq && opt.hastime
    val = opt.fun(xi, yi, zi, opt.qi);
  elseif opt.hasfreq && ~opt.hastime
    val = opt.fun(xi, yi, zi, opt.qi);
  elseif opt.hasfreq && opt.hastime
    val = opt.fun(xi, yi, zi, opt.qi(1), opt.qi(2));
  end
  str4 = sprintf('value %f', val);
else
  str4 = '';
end

%fprintf('%s %s %s %s\n', str1, str2, str3, str4);

if opt.hasatlas
  %tmp = [opt.ijk(:)' 1] * opt.atlas.transform; % atlas and functional might have different transformation matrices, so xyz cannot be used here anymore
  % determine the anatomical label of the current position
  lab = atlas_lookup(opt.atlas, (xyz(1:3)), 'inputcoord', functional.coordsys, 'queryrange', opt.queryrange);
  if isempty(lab)
    lab = 'NA';
    %fprintf('atlas labels: not found\n');
  else
    tmp = sprintf('%s', strrep(lab{1}, '_', ' '));
    for i=2:length(lab)
      tmp = [tmp sprintf(', %s', strrep(lab{i}, '_', ' '))];
    end
    lab = tmp;
  end
else
  lab = 'NA';
end


if opt.hasana
  if opt.init
    tmph  = [h1 h2 h3];
    ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'parents', tmph, 'update', opt.update, 'doscale', false, 'clim', opt.clim);

    opt.anahandles = findobj(opt.handlesfigure, 'type', 'surface')';
    for i=1:length(opt.anahandles)
      opt.parenttag{i} = get(get(opt.anahandles(i), 'parent'), 'tag');
    end
    [i1,i2,i3] = intersect(opt.parenttag, {'ik' 'jk' 'ij'});
    opt.anahandles = opt.anahandles(i3(i2)); % seems like swapping the order
    opt.anahandles = opt.anahandles(:)';
    set(opt.anahandles, 'tag', 'ana');
  else
    ft_plot_ortho(opt.ana, 'transform', eye(4), 'location', opt.ijk, 'style', 'subplot', 'surfhandle', opt.anahandles, 'update', opt.update, 'doscale', false, 'clim', opt.clim);
  end
end

if opt.hasfun
  if opt.init
    if opt.hasmsk
      tmpqi = [opt.qi 1];
      tmph  = [h1 h2 h3];
      ft_plot_ortho(opt.fun(:,:,:,tmpqi(1),tmpqi(2)), 'datmask', opt.msk(:,:,:,tmpqi(1),tmpqi(2)), 'transform', eye(4), 'location', opt.ijk, ...
        'style', 'subplot', 'parents', tmph, 'update', opt.update, ...
        'colormap', opt.funcolormap, 'clim', [opt.fcolmin opt.fcolmax], ...
        'opacitylim', [opt.opacmin opt.opacmax]);


    else
      tmpqi = [opt.qi 1];
      tmph  = [h1 h2 h3];
      ft_plot_ortho(opt.fun(:,:,:,tmpqi(1),tmpqi(2)), 'transform', eye(4), 'location', opt.ijk, ...
        'style', 'subplot', 'parents', tmph, 'update', opt.update, ...
        'colormap', opt.funcolormap, 'clim', [opt.fcolmin opt.fcolmax]);
    end
    % after the first call, the handles to the functional surfaces
    % exist. create a variable containing this, and sort according to
    % the parents
    opt.funhandles = findobj(opt.handlesfigure, 'type', 'surface');
    opt.funtag     = get(opt.funhandles, 'tag');
    opt.funhandles = opt.funhandles(~strcmp('ana', opt.funtag));
    for i=1:length(opt.funhandles)
      opt.parenttag{i} = get(get(opt.funhandles(i), 'parent'), 'tag');
    end
    [i1,i2,i3] = intersect(opt.parenttag, {'ik' 'jk' 'ij'});
    opt.funhandles = opt.funhandles(i3(i2)); % seems like swapping the order
    opt.funhandles = opt.funhandles(:)';
    set(opt.funhandles, 'tag', 'fun');

    if ~opt.hasmsk && opt.hasfun && opt.hasana
      set(opt.funhandles(1), 'facealpha',0.5);
      set(opt.funhandles(2), 'facealpha',0.5);
      set(opt.funhandles(3), 'facealpha',0.5);
    end

  else
    if opt.hasmsk
      tmpqi = [opt.qi 1];
      tmph  = opt.funhandles;
      ft_plot_ortho(opt.fun(:,:,:,tmpqi(1),tmpqi(2)), 'datmask', opt.msk(:,:,:,tmpqi(1),tmpqi(2)), 'transform', eye(4), 'location', opt.ijk, ...
        'style', 'subplot', 'surfhandle', tmph, 'update', opt.update, ...
        'colormap', opt.funcolormap, 'clim', [opt.fcolmin opt.fcolmax], ...
        'opacitylim', [opt.opacmin opt.opacmax]);
    else
      tmpqi = [opt.qi 1];
      tmph  = opt.funhandles;
      ft_plot_ortho(opt.fun(:,:,:,tmpqi(1),tmpqi(2)), 'transform', eye(4), 'location', opt.ijk, ...
        'style', 'subplot', 'surfhandle', tmph, 'update', opt.update, ...
        'colormap', opt.funcolormap, 'clim', [opt.fcolmin opt.fcolmax]);
    end
  end
end
set(opt.handlesaxes(1), 'Visible',opt.axis);
set(opt.handlesaxes(2), 'Visible',opt.axis);
set(opt.handlesaxes(3), 'Visible',opt.axis);

if opt.hasfreq && opt.hastime && opt.hasfun,
  h4 = subplot(2,2,4);
  tmpdat = double(squeeze(opt.fun(xi,yi,zi,:,:)));
  uimagesc(double(functional.time), double(functional.freq), tmpdat); axis xy;
  xlabel('time'); ylabel('freq');
  set(h4, 'tag', 'TF1');
  caxis([opt.fcolmin opt.fcolmax]);
elseif opt.hasfreq && opt.hasfun,
  h4 = subplot(2,2,4);
  plot(functional.freq, squeeze(opt.fun(xi,yi,zi,:))); xlabel('freq');
  axis([functional.freq(1) functional.freq(end) opt.fcolmin opt.fcolmax]);
  set(h4, 'tag', 'TF2');
elseif opt.hastime && opt.hasfun,
  h4 = subplot(2,2,4);
  plot(functional.time, squeeze(opt.fun(xi,yi,zi,:))); xlabel('time');
  set(h4, 'tag', 'TF3', 'xlim',functional.time([1 end]), 'ylim',[opt.fcolmin opt.fcolmax], 'layer', 'top');
elseif strcmp(opt.colorbar,  'yes') && ~isfield(opt, 'hc'),
  if opt.hasfun
    % vectorcolorbar = linspace(fscolmin, fcolmax,length(cfg.funcolormap));
    % imagesc(vectorcolorbar,1,vectorcolorbar);colormap(cfg.funcolormap);
    % use a normal MATLAB colorbar, attach it to the invisible 4th subplot
    try
      caxis([opt.fcolmin opt.fcolmax]);
    end
    opt.hc = colorbar;
    set(opt.hc, 'location', 'southoutside');
    set(opt.hc, 'position',[0.06+0.06+opt.h1size(1) 0.06-0.06+opt.h3size(2) opt.h2size(1) 0.06]);

    try
      set(opt.hc, 'XLim', [opt.fcolmin opt.fcolmax]);
    end
  else
    ft_warning('no colorbar possible without functional data');
  end
end

if ~((opt.hasfreq && numel(functional.freq)>1) || opt.hastime)
  if opt.init
    ht = subplot('position',[0.06+0.06+opt.h1size(1) 0.06 opt.h2size(1) opt.h3size(2)]);
    set(ht, 'visible', 'off');
    opt.ht1=text(0,0.6,str1);
    opt.ht2=text(0,0.5,str2);
    opt.ht3=text(0,0.4,str4);
    opt.ht4=text(0,0.3,str3);
    opt.ht5=text(0,0.2,['atlas label: ' lab]);
  else
    set(opt.ht1, 'string',str1);
    set(opt.ht2, 'string',str2);
    set(opt.ht3, 'string',str4);
    set(opt.ht4, 'string',str3);
    set(opt.ht5, 'string',['atlas label: ' lab]);
  end
end

% make the last current axes current again
sel = findobj('type', 'axes', 'tag',tag);
if ~isempty(sel)
  set(opt.handlesfigure, 'currentaxes', sel(1));
end
if opt.crosshair
  if opt.init
    hch1 = crosshair([xi 1 zi], 'parent', opt.handlesaxes(1));
    hch3 = crosshair([xi yi opt.dim(3)], 'parent', opt.handlesaxes(3));
    hch2 = crosshair([opt.dim(1) yi zi], 'parent', opt.handlesaxes(2));
    opt.handlescross  = [hch1(:)';hch2(:)';hch3(:)'];
  else
    crosshair([xi 1 zi], 'handle', opt.handlescross(1, :));
    crosshair([opt.dim(1) yi zi], 'handle', opt.handlescross(2, :));
    crosshair([xi yi opt.dim(3)], 'handle', opt.handlescross(3, :));
  end
end

if opt.init
  opt.init = false;
  setappdata(h, 'opt', opt);
end

set(h, 'currentaxes', curr_ax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_keyboard(h, eventdata)

if isempty(eventdata)
  % determine the key that corresponds to the uicontrol element that was activated
  key = get(h, 'userdata');
else
  % determine the key that was pressed on the keyboard
  key = parseKeyboardEvent(eventdata);
end
% get focus back to figure
if ~strcmp(get(h, 'type'), 'figure')
  set(h, 'enable', 'off');
  drawnow;
  set(h, 'enable', 'on');
end

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h, 'currentaxes');
tag     = get(curr_ax, 'tag');

if isempty(key)
  % this happens if you press the apple key
  key = '';
end

% the following code is largely shared with FT_VOLUMEREALIGN
switch key
  case {'' 'shift+shift' 'alt-alt' 'control+control' 'command-0'}
    % do nothing

  case '1'
    subplot(opt.handlesaxes(1));

  case '2'
    subplot(opt.handlesaxes(2));

  case '3'
    subplot(opt.handlesaxes(3));

  case 'q'
    setappdata(h, 'opt', opt);
    cb_cleanup(h);

  case {'i' 'j' 'k' 'm' 28 29 30 31 'leftarrow' 'rightarrow' 'uparrow' 'downarrow'} % TODO FIXME use leftarrow rightarrow uparrow downarrow
    % update the view to a new position
    if     strcmp(tag, 'ik') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.ijk(3) = opt.ijk(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag, 'ik') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.ijk(1) = opt.ijk(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ik') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.ijk(1) = opt.ijk(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ik') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.ijk(3) = opt.ijk(3)-1; opt.update = [0 0 1];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.ijk(2) = opt.ijk(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.ijk(1) = opt.ijk(1)-1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.ijk(1) = opt.ijk(1)+1; opt.update = [0 1 0];
    elseif strcmp(tag, 'ij') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.ijk(2) = opt.ijk(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'i') || strcmp(key, 'uparrow')    || isequal(key, 30)), opt.ijk(3) = opt.ijk(3)+1; opt.update = [0 0 1];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'j') || strcmp(key, 'leftarrow')  || isequal(key, 28)), opt.ijk(2) = opt.ijk(2)-1; opt.update = [1 0 0];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'k') || strcmp(key, 'rightarrow') || isequal(key, 29)), opt.ijk(2) = opt.ijk(2)+1; opt.update = [1 0 0];
    elseif strcmp(tag, 'jk') && (strcmp(key, 'm') || strcmp(key, 'downarrow')  || isequal(key, 31)), opt.ijk(3) = opt.ijk(3)-1; opt.update = [0 0 1];
    else
      % do nothing
    end;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

    % contrast scaling
  case {43 'shift+equal'}  % numpad +
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % reduce color scale range by 5%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
    %opt.clim(1) = opt.clim(1)+cscalefactor;
    opt.clim(2) = opt.clim(2)-cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

  case {45 'shift+hyphen'} % numpad -
    if isempty(opt.clim)
      opt.clim = [min(opt.ana(:)) max(opt.ana(:))];
    end
    % increase color scale range by 5%
    cscalefactor = (opt.clim(2)-opt.clim(1))/10;
    %opt.clim(1) = opt.clim(1)-cscalefactor;
    opt.clim(2) = opt.clim(2)+cscalefactor;
    setappdata(h, 'opt', opt);
    cb_redraw(h);

  otherwise

end % switch key

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonpress(h, eventdata)

h   = getparent(h);
cb_getposition(h);

switch get(h, 'selectiontype')
  case 'normal'
    % just update to new position, nothing else to be done here
    cb_redraw(h);
  case 'alt'
    set(h, 'windowbuttonmotionfcn', @cb_tracemouse);
    cb_redraw(h);
  otherwise
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_buttonrelease(h, eventdata)

set(h, 'windowbuttonmotionfcn', '');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_tracemouse(h, eventdata)

h   = getparent(h);
cb_getposition(h);
cb_redraw(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_getposition(h, eventdata)

h   = getparent(h);
opt = getappdata(h, 'opt');

curr_ax = get(h,       'currentaxes');
pos     = mean(get(curr_ax, 'currentpoint'));

tag = get(curr_ax, 'tag');

if ~isempty(tag) && ~opt.init
  if strcmp(tag, 'ik')
    opt.ijk([1 3])  = round(pos([1 3]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'ij')
    opt.ijk([1 2])  = round(pos([1 2]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'jk')
    opt.ijk([2 3])  = round(pos([2 3]));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'TF1')
    % timefreq
    opt.qi(2) = nearest(opt.functional.time, pos(1));
    opt.qi(1) = nearest(opt.functional.freq, pos(2));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'TF2')
    % freq only
    opt.qi  = nearest(opt.functional.freq, pos(1));
    opt.update = [1 1 1];
  elseif strcmp(tag, 'TF3')
    % time only
    opt.qi  = nearest(opt.functional.time, pos(1));
    opt.update = [1 1 1];
  end
end
opt.ijk = min(opt.ijk(:)', opt.dim);
opt.ijk = max(opt.ijk(:)', [1 1 1]);

setappdata(h, 'opt', opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cb_cleanup(h, eventdata)

% opt = getappdata(h, 'opt');
% opt.quit = true;
% setappdata(h, 'opt', opt);
% uiresume
delete(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = getparent(h)
p = h;
while p~=0
  h = p;
  p = get(h, 'parent');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = parseKeyboardEvent(eventdata)

key = eventdata.Key;

% handle possible numpad events (different for Windows and UNIX systems)
% NOTE: shift+numpad number does not work on UNIX, since the shift
% modifier is always sent for numpad events
if isunix()
  shiftInd = match_str(eventdata.Modifier, 'shift');
  if ~isnan(str2double(eventdata.Character)) && ~isempty(shiftInd)
    % now we now it was a numpad keystroke (numeric character sent AND
    % shift modifier present)
    key = eventdata.Character;
    eventdata.Modifier(shiftInd) = []; % strip the shift modifier
  end
elseif ispc()
  if strfind(eventdata.Key, 'numpad')
    key = eventdata.Character;
  end
end

if ~isempty(eventdata.Modifier)
  key = [eventdata.Modifier{1} '+' key];
end
