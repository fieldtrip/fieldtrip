function [data] = ft_megplanar(cfg, data)

% FT_MEGPLANAR computes planar MEG gradients gradients for raw data or average
% event-related field data. It can also convert frequency-domain data that was computed
% using FT_FREQANALYSIS, as long as it contains the complex-valued fourierspcrm and not
% only the powspctrm.
%
% Use as
%    [interp] = ft_megplanar(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING,
% FT_TIMELOCKANALYSIS or FT_FREQANALYSIS (with output='fourierspcrm').
%
% The configuration should contain
%   cfg.planarmethod   = string, can be 'sincos', 'orig', 'fitplane', 'sourceproject' (default = 'sincos')
%   cfg.channel        =  Nx1 cell-array with selection of channels (default = 'MEG'), see FT_CHANNELSELECTION for details
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%
% The methods orig, sincos and fitplane are all based on a neighbourhood interpolation.
% For these methods you need to specify
%   cfg.neighbours     = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%
% In the 'sourceproject' method a minumum current estimate is done using a large number
% of dipoles that are placed in the upper layer of the brain surface, followed by a
% forward computation towards a planar gradiometer array. This requires the
% specification of a volume conduction model of the head and of a source model. The
% 'sourceproject' method is not supported for frequency domain data.
%
% A dipole layer representing the brain surface must be specified with
%   cfg.inwardshift = depth of the source layer relative to the head model surface (default = 2.5 cm, which is appropriate for a skin-based head model)
%   cfg.spheremesh  = number of dipoles in the source layer (default = 642)
%   cfg.pruneratio  = for singular values, default is 1e-3
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
% If no headshape is specified, the dipole layer will be based on the inner compartment
% of the volume conduction model.
%
% The volume conduction model of the head should be specified as
%   cfg.vol           = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%   cfg.hdmfile       = name of file containing the volume conduction model, see FT_READ_VOL
%
% The following cfg fields are optional:
%   cfg.feedback
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
% See also FT_COMBINEPLANAR, FT_NEIGHBOURSELECTION

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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble provenance
ft_preamble trackconfig
ft_preamble loadvar data

% store the original input representation of the data, this is used later on to convert it back
isfreq = ft_datatype(data, 'freq');
israw  = ft_datatype(data, 'raw');
istlck = ft_datatype(data, 'timelock');  % this will be temporary converted into raw

% check if the input data is valid for this function, this converts the data if needed
data = ft_checkdata(data, 'datatype', {'raw' 'freq'}, 'feedback', 'yes', 'hassampleinfo', 'yes', 'ismeg', 'yes', 'senstype', {'ctf151', 'ctf275', 'bti148', 'bti248', 'itab153', 'yokogawa160', 'yokogawa64'});

if istlck
  % the timelocked data has just been converted to a raw representation
  % and will be converted back to timelocked at the end of this function
  israw = true;
end

if isfreq,
  if ~isfield(data, 'fourierspctrm'), error('freq data should contain Fourier spectra'); end
end

% set the default configuration
cfg.channel      = ft_getopt(cfg, 'channel',      'MEG');
cfg.trials       = ft_getopt(cfg, 'trials',       'all');
cfg.planarmethod = ft_getopt(cfg, 'planarmethod', 'sincos');
cfg.inputfile    = ft_getopt(cfg, 'inputfile',    []);
cfg.outputfile   = ft_getopt(cfg, 'outputfile',   []);
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');

% check if the input cfg is valid for this function
if ~strcmp(cfg.planarmethod, 'sourceproject')
  cfg = ft_checkconfig(cfg, 'required', {'neighbours'});
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

if strcmp(cfg.planarmethod, 'sourceproject')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do an inverse computation with a simplified distributed source model
  % and compute forward again with the axial gradiometer array replaced by
  % a planar one.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % method specific configuration options
  cfg.headshape   = ft_getopt(cfg, 'headshape',   []);
  cfg.inwardshift = ft_getopt(cfg, 'inwardshift', 2.5); % this number assumes that all other inputs are in cm
  cfg.pruneratio  = ft_getopt(cfg, 'pruneratio',  1e-3);
  cfg.spheremesh  = ft_getopt(cfg, 'spheremesh',  642);
  
  if isfreq
    error('the method ''sourceproject'' is not supported for frequency data as input');
  end
  
  Nchan   = length(data.label);
  Ntrials = length(data.trial);
  
  % FT_PREPARE_VOL_SENS will match the data labels, the gradiometer labels and the
  % volume model labels (in case of a localspheres model) and result in a gradiometer
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
  % this will not work for a localspheres model, compute_leadfield will catch
  % the error
  planar.grad = constructplanargrad([], axial.grad);
  lfnew = ft_compute_leadfield(pos, planar.grad, vol);
  
  % compute the interpolation matrix
  transform = lfnew * prunedinv(lfold, cfg.pruneratio);
  
  planarmontage = [];
  planarmontage.tra      = transform;
  planarmontage.labelorg = axial.grad.label;
  planarmontage.labelnew = planar.grad.label;
  
  % apply the linear transformation to the data
  interp  = ft_apply_montage(data, planarmontage, 'keepunused', 'yes');
  % also apply the linear transformation to the gradiometer definition
  interp.grad = ft_apply_montage(data.grad, planarmontage, 'balancename', 'planar', 'keepunused', 'yes');
  
  % ensure there is a type string describing the gradiometer definition
  if ~isfield(interp.grad, 'type')
    % put the original gradiometer type in (will get _planar appended)
    interp.grad.type = ft_senstype(data.grad);
  end
  interp.grad.type = [interp.grad.type '_planar'];
  
  %   % interpolate the data towards the planar gradiometers
  %   for i=1:Ntrials
  %     fprintf('interpolating trial %d to planar gradiometer\n', i);
  %     interp.trial{i} = transform * data.trial{i}(dataindx,:);
  %   end % for Ntrials
  %
  %   % all planar gradiometer channels are included in the output
  %   interp.grad  = planar.grad;
  %   interp.label = planar.grad.label;
  %
  %   % copy the non-gradiometer channels back into the output data
  %   other = setdiff(1:Nchan, dataindx);
  %   for i=other
  %     interp.label{end+1} = data.label{i};
  %     for j=1:Ntrials
  %       interp.trial{j}(end+1,:) = data.trial{j}(i,:);
  %     end
  %   end
  %
else
  % generically call megplanar_orig megplanar_sincos or megplanar_fitplante
  fun = ['megplanar_'  cfg.planarmethod];
  if ~exist(fun, 'file')
    error('unknown method for computation of planar gradient');
  end
  
  sens = ft_convert_units(data.grad);
  if any(isnan(sens.chanpos(:)))
    error('The channel positions contain NaNs; this prohibits correct behavior of the function. Please replace the input channel definition with one that contains valid channel positions');
  end
  cfg.channel = ft_channelselection(cfg.channel, sens.label);
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  
  % ensure channel order according to cfg.channel (there might be one check
  % too much in here somewhere or in the subfunctions, but I don't care.
  % Better one too much than one too little - JMH @ 09/19/12
  [neighbsel] = match_str({cfg.neighbours.label}, cfg.channel);
  cfg.neighbours = cfg.neighbours(neighbsel);
  cfg.neighbsel = channelconnectivity(cfg);
  
  % determine
  fprintf('average number of neighbours is %.2f\n', mean(sum(cfg.neighbsel)));
  
  Ngrad = length(sens.label);
  cfg.distance = zeros(Ngrad,Ngrad);
  
  for i=1:size(cfg.neighbsel,1)
    j=find(cfg.neighbsel(i, :));
    d = sqrt(sum((sens.chanpos(j,:) - repmat(sens.chanpos(i, :), numel(j), 1)).^2, 2));
    distance(i,j) = d;
    distance(j,i) = d;
  end
  
  fprintf('minimum distance between neighbours is %6.2f %s\n', min(distance(distance~=0)), sens.unit);
  fprintf('maximum distance between gradiometers is %6.2f %s\n', max(distance(distance~=0)), sens.unit);
    
  planarmontage = eval([fun '(cfg, data.grad)']);
  
  % apply the linear transformation to the data
  interp = ft_apply_montage(data, planarmontage, 'keepunused', 'yes', 'feedback', cfg.feedback);
  
  % also apply the linear transformation to the gradiometer definition
  interp.grad = ft_apply_montage(data.grad, planarmontage, 'balancename', 'planar', 'keepunused', 'yes');
  
  % ensure there is a type string describing the gradiometer definition
  if ~isfield(interp.grad, 'type')
    % put the original gradiometer type in (will get _planar appended)
    interp.grad.type = ft_senstype(data.grad);
  end
  interp.grad.type = [interp.grad.type '_planar'];
  
  % add the chanpos info back into the gradiometer description
  tmplabel = interp.grad.label;
  for k = 1:numel(tmplabel)
    if strcmp(tmplabel{k}(end-2:end), '_dV') || strcmp(tmplabel{k}(end-2:end), '_dH')
      tmplabel{k} = tmplabel{k}(1:end-3);
    end
  end
  [ix,iy] = match_str(tmplabel, data.grad.label);
  interp.grad.chanpos(ix,:) = data.grad.chanpos(iy,:);
end

if istlck
  % convert the raw structure back into a timelock structure
  interp = ft_checkdata(interp, 'datatype', 'timelock');
  israw  = false;
end

% copy the trial specific information into the output
if isfield(data, 'trialinfo')
  interp.trialinfo = data.trialinfo;
end

% copy the sampleinfo field as well
if isfield(data, 'sampleinfo')
  interp.sampleinfo = data.sampleinfo;
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble provenance
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = interp;

ft_postamble history data
ft_postamble savevar data


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
