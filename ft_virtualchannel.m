function [data_vc] = ft_virtualchannel(cfg, data, source, parcellation)

% FT_VIRTUALCHANNEL creates virtual channel data, combining numeric data 
% from a data structure defined at the channel level with spatial filter
% information from a source data structure, and optional parcellation
% information.
%
% Use as
%    output = ft_virtualchannel(cfg, data, source)
% or 
%    output = ft_virtualchannel(cfg, data, source, parcellation)
%
% where the input "data" is a channel level data that contains data that can
% be linearly mapped onto the virtual channel level, e.g. a raw data
% structure obtained with FT_PREPROCESSING, a timelock structure, obtained
% with FT_TIMELOCKANALYSIS, or a freq structure with fourierspectra,
% obtained with FT_FREQANALYSIS. The input "source" is a source structure
% that has been obtained with FT_SOURCEANALYSIS, and which contains spatial
% filter information for at least one dipole location, in the
% source.filter, or source.avg.filter field. The optional input
% "parcellation" is described in detail in FT_DATATYPE_PARCELLATION (2-D) or 
% FT_DATATYPE_SEGMENTATION (3-D) and can be obtained from FT_READ_ATLAS or
% from a custom parcellation/segmentation for your individual subject.
% Alternatively, the input "source" can already contain a parcellation.
%
% The configuration "cfg" is a structure that should either contain
%   cfg.pos    = Nx3 matrix containing the dipole positions for the virtual
%                  channel(s). These positions should match the entries in
%                  the source.pos field. (default = [])
% or
%   cfg.parcellation = string, name of the field that is used for the
%                  parcel labels. (default = [])
%   cfg.parcel = string, or cell-array of strings, specifying for which
%                  parcels to return the output. (default = 'all')
%
% The values within a parcel or parcel-combination can be combined with different methods:
%   'mean'      compute the mean
%   'median'    compute the median (unsupported for fields that are represented in a cell-array)
%   'eig'       compute the largest eigenvector
%   'min'       take the minimal value
%   'max'       take the maximal value
%   'maxabs'    take the signed maxabs value
%   'std'       take the standard deviation
%
% See also FT_SOURCEANALYSIS, FT_DATATYPE_PARCELLATION, FT_DATATYPE_SEGMENTATION

% Copyright (C) 2020, Jan-Mathijs Schoffelen and Robert Oostenveld
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
ft_preamble loadvar data source parcellation
ft_preamble provenance data source parcellation
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% get the defaults
cfg.pos          = ft_getopt(cfg, 'pos');
cfg.parcellation = ft_getopt(cfg, 'parcellation');
cfg.parcel       = ft_getopt(cfg, 'parcel',   'all');
cfg.method       = ft_getopt(cfg, 'method',   'svd');
cfg.feedback     = ft_getopt(cfg, 'feedback', 'text');
cfg.numcomponent = ft_getopt(cfg, 'numcomponent', 1);

% the data can be passed as input argument or can be read from disk
hasparcellation = exist('parcellation', 'var');

if hasparcellation
  % the parcellation is specified as separate structure
elseif ~isempty(cfg.parcellation)
  % the parcellation is represented in the source structure itself
  parcellation = source;
end

useparcellation = ~isempty(cfg.parcellation) && isfield(parcellation, cfg.parcellation);
usepos          = ~isempty(cfg.pos);

if usepos && useparcellation
  ft_error('you should either specify cfg.pos, or cfg.parcellation, not both');
elseif usepos
  % extract the virtual channels from the specified positions
elseif useparcellation
  % extract the virtual channels from the specified parcels
else
  ft_error('you should either specify cfg.pos, or a valid cfg.parcellation');
end

% ensure that the source input is a source, not a volume, and this should also return the
% source.filter, rather than source.avg.filter
source = ft_checkdata(source, 'datatype', 'source', 'inside', 'logical', 'hasunit', 'yes');
if ~isfield(source, 'filter')
  ft_error('the input source needs a ''filter'' field');
  % FIXME how about the output of ft_dipolefitting?
end

% FIXME for now support only raw-type data structures, this should be
% extended to also allow freq (with fourier) and timelock structures,
% although the latter is probably already covered by the call to checkdata,
% it should be checked how this is done in ft_megplanar.
data     = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'mvar'}, 'feedback', cfg.feedback);
[i1, i2] = match_str(data.label, source.label);

% select channels and trials of interest, by default this will select all channels and trials
tmpcfg = keepfields(cfg, {'trials', 'tolerance', 'showcallinfo'});
tmpcfg.channel = data.label(i1);
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% do the actual work
if usepos
  
  % identify the indices of the positions in the source structure that are
  % to be used, and check whether a spatial filter exists
  npos  = size(source.pos,1);
  nvc   = size(cfg.pos,1);
  indx  = zeros(nvc,1);
  mindist = zeros(nvc,1);
  for i = 1:nvc
    dpos = sqrt(sum( (source.pos - cfg.pos(i.*ones(npos,1),:)).^2, 2));
    [mindist(i), indx(i)] = min(dpos);
    
    % check that the requested positions are at most 1 mm away from the actual
    % positions
    if mindist(i)>1*ft_scalingfactor('mm', source.unit)
      ft_error('the requested dipole position is > 1 mm away from the closest source position');
    end
    
    % check that the spatial filter exists
    if isempty(source.filter{indx(i)})
      ft_error(sprintf('the spatial filter information is missing for dipole with position [%d %d %d]', cfg.pos(i,:)));
    end
  end
  
elseif useparcellation
  
  % do some checks and balances on the input data and cfg
  
  % keep the transformation matrix
  if isfield(parcellation, 'transform')
    transform = parcellation.transform;
  else
    transform = [];
  end
  
  % ensure it is a parcellation, not a segmentation
  parcellation = ft_checkdata(parcellation, 'datatype', 'parcellation', 'parcel lationstyle', 'indexed');
  
  % keep the transformation matrix
  if ~isempty(transform)
    parcellation.transform = transform;
  end
 
  % ensure that the source and the parcellation are anatomically consistent
  if ~isalmostequal(source.pos, parcellation.pos, 'abstol', 1000000*eps)
    ft_error('the source positions are not consistent with the parcellation, please use FT_SOURCEINTERPOLATE');
  end
  
  if isequal(cfg.parcel, 'all')
    cfg.parcel = parcellation.(sprintf('%slabel',cfg.parcellation));
  end
  
  indx = match_str(parcellation.(sprintf('%slabel',cfg.parcellation)), cfg.parcel);
  cfg.parcel = parcellation.(sprintf('%slabel',cfg.parcellation))(indx);
  
  nvc = numel(cfg.parcel);
end

% create a montage for each of the dipoles/parcels, and use ft_apply_montage,
% followed by a dimensionality reduction step, using ft_componentanalysis
unmixing  = cell(1, nvc);
topolabel = cell(1, nvc);
label     = cell(1, nvc);

compmethods = {'pca' 'runica' 'fastica' 'dss'};
for i = 1:nvc
  
  montage = [];
  if usepos
    
    montage.tra      = source.filter{indx(i)};
    montage.labelold = source.label;
    montage.labelnew = cell(size(montage.tra,1), 1);
    for k = 1:size(montage.tra,1)
      montage.labelnew{k} = sprintf('virtualchannel%03d_orientation%03d', i, k);
    end
    bname = sprintf('virtualchannel%03d', i);
    
  elseif useparcellation
    
    sel = find(parcellation.(cfg.parcellation)==indx(i));
    montage.tra      = cat(1, source.filter{sel});
    montage.labelold = source.label;
    montage.labelnew = cell(size(montage.tra,1), 1);
    cnt = 0;
    for k = 1:numel(sel)
      for kk = 1:size(source.filter{sel(k)},1)
        cnt = cnt+1;
        montage.labelnew{cnt} = sprintf('%s_dipole%03d_orientation%03d', cfg.parcel{i}, k, kk);
      end
    end
    bname = matlab.lang.makeValidName(cfg.parcel{i}); % FIXME not sure how universal this is.
    
  end
  
  [i1, i2] = match_str(data.label, montage.labelold);
  
  % ensure that the used channels run from 1:numel(data.label), to be
  % sure that the to-be-computed covariance matches the labels, not sure
  % whether this is needed, this could probably also be achieved by doing
  % a check on the input data's labels, related to the source.filter
  % labels
  assert(isequal(i1(:)', 1:numel(data.label)));
  
  unmixing{1,i}(:,i1) = montage.tra(:,i2);
  topolabel{1,i}      = data.label(i1);
  
  % Next, apply the montage to the numeric data, but make an exception for the
  % svd method: this one is much more efficient, in that its implementation
  % operates on the covariance matrix directly, bypasses ft_componentanalysis,
  % and does not require the intermediate (possibly memory greedy, and slow)
  % data representation.
  switch cfg.method
    case 'svd'
          
      if ~exist('tlck', 'var')
        % create a data structure that contains a covariance -> FIXME requires
        % a check on the input data, for now assume raw data
        tmpcfg = [];
        tmpcfg.covariance = 'yes';
        tlck   = ft_timelockanalysis(tmpcfg, data);
      end
      
      % create a matricial square root of the full covariance matrix
      [u, s, v] = svd(tlck.cov);
      Csqrtm    = u*sqrt(s);
      
      % don't do the sandwiching for efficiency, we will use only the u
      % matrix, this yields the same u matrix as the svd on the tlck.cov
      [u, s, v] = svd(unmixing{1, i} * Csqrtm, 'econ');
      
      if isequal(cfg.numcomponent, 'all')
        ncomp = size(u,2);
      else
        if size(u,2)>=cfg.numcomponent
          ncomp = cfg.numcomponent;
        else
          ncomp = size(u,2);
          ft_warning(sprintf('using %d components,rather than the requested %d\n', ncomp, cfg.numcomponent));
        end
      end
      
      unmixing{1, i} = u(:,1:ncomp)' * unmixing{1, i};
      label{1, i}    = cell(ncomp,1);
      if usepos
        str = sprintf('virtualchannel%03d', i);
      else
        str = cfg.parcel{i};
      end
      if ncomp == 1
        label{1, i}{1} = str;
      else
        for k = 1:ncomp
          label{1, i}{k} = sprintf('%s_%03d', str, k);
        end
      end
      
    case 'none'
      % do nothing, i.e. stick to the original 'montage', which is the
      % 3D, or concatenated spatial filter
      
    case compmethods
      
      % apply the montage to the data
      tmpdata = ft_apply_montage(data, montage, 'feedback', 'none');
      
      tmpcfg     = keepfields(cfg, {'method' cfg.method 'numcomponent'});
      if isequal(tmpcfg.numcomponent, 'all')
        tmpcfg = rmfield(tmpcfg, 'numcomponent');
      end
      
      tmpdata     = ft_componentanalysis(tmpcfg, tmpdata);
      unmixing{i} = tmpdata.unmixing * unmixing{i};
      
      for k = 1:numel(tmpdata.label)
        if usepos
          str = sprintf('virtualchannel%03d', i);
        else
          str = cfg.channel{i};
        end
        if numel(tmpdata.label)==1
          tmpdata.label{k} = sprintf('%s_%03d', str);
        else
          tmpdata.label{k} = sprintf('%s_%03d', str, k);
        end
      end
      
    otherwise
      ft_error('currently not yet supported');
      % the idea would to support a custom function here, with a
      % function(cfg.(cfg.method), tmpdata{i}) API
      
  end % reduction of components
  
end % for i = # of virtual channels

% create and apply the compound montage
montage          = [];
montage.tra      = cat(1, unmixing{:});
montage.labelnew = cat(1, label{:});
montage.labelold = topolabel{1};

data_vc = ft_apply_montage(data, montage, 'feedback', 'none');

% apply the montage to the sensor description
sensfields = {'grad' 'elec' 'opto'};
bname      = cfg.parcellation;
for k = 1:numel(sensfields)
  if isfield(data_vc, sensfields{k})
    ft_info(sprintf('applying the montage to the %s structure\n', sensfields{k}));
    data_vc.(sensfields{k}) = ft_apply_montage(data.grad, montage, 'feedback', 'none', 'keepunused', 'no', 'balancename', bname);
  end
end

if usepos
  brainordinate = keepfields(source, {'pos' 'dim' 'tri' 'transform' 'inside' 'unit'});
  brainordinate.index = zeros(npos, 1);
  brainordinate.index(indx) = indx;
  brainordinate.indexlabel  = data_vc.label; % FIXME this only works with one component per vc
elseif useparcellation
  brainordinate = parcellation;
end
data_vc.brainordinate = brainordinate;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data source parcellation 
ft_postamble provenance data_vc
ft_postamble history    data_vc
ft_postamble savevar    data_vc
