function [interp] = ft_megplanar(cfg, data)

% FT_MEGPLANAR computes planar MEG gradients gradients for raw data
% obtained from FT_PREPROCESSING or an average ERF that was computed using
% FT_TIMELOCKANALYSIS. It can also work on data in the frequency domain,
% obtained with FREQANALYSIS. Prerequisite for this is that the data contain
% complex-valued fourierspectra.
%
% Use as
%    [interp] = ft_megplanar(cfg, data)
%
% The configuration should contain
%   cfg.planarmethod   = 'orig' | 'sincos' | 'fitplane' | 'sourceproject'
%   cfg.channel        =  Nx1 cell-array with selection of channels (default = 'MEG'),
%                         see FT_CHANNELSELECTION for details
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%
% The methods orig, sincos and fitplane are all based on a neighbourhood
% interpolation. For these methods you can specify
%   cfg.neighbourdist  = default is 4 cm
%
% In the 'sourceproject' method a minumum current estimate is done using a
% large number of dipoles that are placed in the upper layer of the brain
% surface, followed by a forward computation towards a planar gradiometer
% array. This requires the specification of a volume conduction model of
% the head and of a source model. The 'sourceproject' method is not supported for
% frequency domain data.
%
% A head model must be specified with
%   cfg.hdmfile     = string, file containing the volume conduction model
% or alternatively manually using
%   cfg.vol.r       = radius of sphere
%   cfg.vol.o       = [x, y, z] position of origin
%
% A dipole layer representing the brain surface must be specified with
%   cfg.inwardshift = depth of the source layer relative to the head model surface (default = 2.5, which is adequate for a skin-based head model)
%   cfg.spheremesh  = number of dipoles in the source layer (default = 642)
%   cfg.pruneratio  = for singular values, default is 1e-3
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
% If no headshape is specified, the dipole layer will be based on the inner compartment
% of the volume conduction model.
%
% See also FT_COMBINEPLANAR

% Undocumented local options:
%   cfg.inputfile        = one can specifiy preanalysed saved data as input
%   cfg.outputfile       = one can specify output as file to save to disk

% This function depends on FT_PREPARE_BRAIN_SURFACE which has the following options:
% cfg.headshape  (default set in FT_MEGPLANAR: cfg.headshape = 'headmodel'), documented
% cfg.inwardshift (default set in FT_MEGPLANAR: cfg.inwardshift = 2.5), documented
% cfg.spheremesh (default set in FT_MEGPLANAR: cfg.spheremesh = 642), documented
%
% This function depends on FT_PREPARE_VOL_SENS which has the following options:
% cfg.channel
% cfg.elec
% cfg.elecfile
% cfg.grad
% cfg.gradfile
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented

% Copyright (C) 2004, Robert Oostenveld
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

ft_defaults

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

% set defaults
if ~isfield(cfg, 'inputfile'),      cfg.inputfile  = [];          end
if ~isfield(cfg, 'outputfile'),     cfg.outputfile = [];          end

% load optional given inputfile as data
hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    data = loadvar(cfg.inputfile, 'data');
  end
end

isfreq = ft_datatype(data, 'freq');
israw  = ft_datatype(data, 'raw');
istlck = ft_datatype(data, 'timelock');  % this will be temporary converted into raw

% check if the input data is valid for this function
data  = ft_checkdata(data, 'datatype', {'raw' 'freq'}, 'feedback', 'yes', 'hastrialdef', 'yes', 'ismeg', 'yes', 'senstype', {'ctf151', 'ctf275', 'bti148', 'bti248', 'itab153'});

if istlck
  % the timelocked data has just been converted to a raw representation
  % and will be converted back to timelocked at the end of this function
  israw = true;
end

if isfreq,
  if ~isfield(data, 'fourierspctrm'), error('freq data should contain Fourier spectra'); end
end

% set the default configuration
if ~isfield(cfg, 'channel'),       cfg.channel = 'MEG';             end
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';              end
% use a smart default for the distance
if ~isfield(cfg, 'neighbourdist'),
  if     isfield(data.grad, 'unit') && strcmp(data.grad.unit, 'cm')
    cfg.neighbourdist = 4;
  elseif isfield(data.grad, 'unit') && strcmp(data.grad.unit, 'mm')
    cfg.neighbourdist = 40;
  else
    % don't provide a default in case the dimensions of the sensor array are unknown
  end
end
if ~isfield(cfg, 'planarmethod'),  cfg.planarmethod = 'sincos';     end
if strcmp(cfg.planarmethod, 'sourceproject')
  if ~isfield(cfg, 'headshape'),     cfg.headshape = [];            end % empty will result in the vol being used
  if ~isfield(cfg, 'inwardshift'),   cfg.inwardshift = 2.5;         end
  if ~isfield(cfg, 'pruneratio'),    cfg.pruneratio = 1e-3;         end
  if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 642;          end
end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'grid'});
cfg = ft_checkconfig(cfg, 'renamedvalue',  {'headshape', 'headmodel', []});

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);
end

if     strcmp(cfg.planarmethod, 'orig')
  montage = megplanar_orig(cfg, data.grad);
  
elseif strcmp(cfg.planarmethod, 'sincos')
  montage = megplanar_sincos(cfg, data.grad);
  
elseif strcmp(cfg.planarmethod, 'fitplane')
  montage = megplanar_fitplane(cfg, data.grad);
  
elseif strcmp(cfg.planarmethod, 'sourceproject')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do an inverse computation with a simplified distributed source model
  % and compute forward again with the axial gradiometer array replaced by
  % a planar one.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if isfreq
    error('the method ''sourceproject'' is not supported for frequency data as input');
  end
  
  Nchan   = length(data.label);
  Ntrials = length(data.trial);
  
  % FT_PREPARE_VOL_SENS will match the data labels, the gradiometer labels and the
  % volume model labels (in case of a multisphere model) and result in a gradiometer
  % definition that only contains the gradiometers that are present in the data.
  [vol, axial.grad, cfg] = prepare_headmodel(cfg, data);
  
  % determine the dipole layer that represents the surface of the brain
  if isempty(cfg.headshape)
    % construct from the inner layer of the volume conduction model
    pos = headsurface(vol, axial.grad, 'surface', 'cortex', 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  else
    % get the surface describing the head shape
    if isstruct(cfg.headshape) && isfield(cfg.headshape, 'pnt')
      % use the headshape surface specified in the configuration
      headshape = cfg.headshape;
    elseif isnumeric(cfg.headshape) && size(cfg.headshape,2)==3
      % use the headshape points specified in the configuration
      headshape.pnt = cfg.headshape;
    elseif ischar(cfg.headshape)
      % read the headshape from file
      headshape = ft_read_headshape(cfg.headshape);
    else
      error('cfg.headshape is not specified correctly')
    end
    if ~isfield(headshape, 'tri')
      % generate a closed triangulation from the surface points
      headshape.pnt = unique(headshape.pnt, 'rows');
      headshape.tri = projecttri(headshape.pnt);
    end
    % construct from the head surface
    pos = headsurface([], [], 'headshape', headshape, 'inwardshift', cfg.inwardshift, 'npnt', cfg.spheremesh);
  end
  
  % compute the forward model for the axial gradiometers
  fprintf('computing forward model for %d dipoles\n', size(pos,1));
  lfold = ft_compute_leadfield(pos, axial.grad, vol);
  
  % construct the planar gradient definition and compute its forward model
  % this will not work for a multisphere model, compute_leadfield will catch
  % the error
  planar.grad = constructplanargrad([], axial.grad);
  lfnew = ft_compute_leadfield(pos, planar.grad, vol);
  
  % compute the interpolation matrix
  transform = lfnew * prunedinv(lfold, cfg.pruneratio);
  
  % interpolate the data towards the planar gradiometers
  for i=1:Ntrials
    fprintf('interpolating trial %d to planar gradiometer\n', i);
    interp.trial{i} = transform * data.trial{i}(dataindx,:);
  end % for Ntrials
  
  % all planar gradiometer channels are included in the output
  interp.grad  = planar.grad;
  interp.label = planar.grad.label;
  
  % copy the non-gradiometer channels back into the output data
  other = setdiff(1:Nchan, dataindx);
  for i=other
    interp.label{end+1} = data.label{i};
    for j=1:Ntrials
      interp.trial{j}(end+1,:) = data.trial{j}(i,:);
    end
  end
  
else
  error('unknown method for computation of planar gradient');
end % cfg.planarmethod

if any(strcmp(cfg.planarmethod, {'orig', 'sincos', 'fitplane'}))
  % apply the linear transformation to the data
  interp  = ft_apply_montage(data, montage, 'keepunused', 'yes');
  % also apply the linear transformation to the gradiometer definition
  interp.grad = ft_apply_montage(data.grad, montage);
  interp.grad.balance.planar = montage;
  % ensure that the old sensor type does not stick around, because it is now invalid
  % the sensor type is added in FT_PREPARE_VOL_SENS but is not used in external fieldtrip code
  if isfield(interp.grad, 'type')
    interp.grad = rmfield(interp.grad, 'type');
  end
  % keep track of the linear transformations that have been applied
  if strcmp(interp.grad.balance.current, 'none')
    interp.grad.balance.current = 'planar';
  else
    interp.grad.balance.current = ['planar_' interp.grad.balance.current];
  end
end

if istlck
  % convert the raw structure back into a timelock structure
  interp = ft_checkdata(interp, 'datatype', 'timelock');
  israw  = false;
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes');

% store the configuration of this function call, including that of the previous function call
cfg.version.name = mfilename('fullpath');
cfg.version.id   = '$Id$';

% add information about the Matlab version used to the configuration
cfg.version.matlab = version();

% remember the configuration details of the input data
try cfg.previous = data.cfg; end

% remember the exact configuration details in the output
interp.cfg = cfg;

% copy the trial specific information into the output
if isfield(data, 'trialinfo')
  interp.trialinfo = data.trialinfo;
end

% copy the sampleinfo field as well
if isfield(data, 'sampleinfo')
  interp.sampleinfo = data.sampleinfo;
end

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', interp); % use the variable name "data" in the output file
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that computes the inverse using a pruned SVD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lfi] = prunedinv(lf, r)
[u, s, v] = svd(lf);
p = find(s<(s(1,1)*r) & s~=0);
fprintf('pruning %d out of %d singular values\n', length(p), min(size(s)));
s(p) = 0;
s(find(s~=0)) = 1./s(find(s~=0));
lfi = v * s' * u';
