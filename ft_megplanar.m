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
%   cfg.inwardshift = depth of the source layer relative to the head model surface ,
%                     (default = 2.5 cm, which is appropriate for a skin-based head model)
%   cfg.spheremesh  = number of dipoles in the source layer (default = 642)
%   cfg.tolerance   = tolerance ratio for leadfield matrix inverse based on a truncated svd,
%                     reflects the relative magnitude of the largest singular value
%                     to retain (default = 1e-3)
%   cfg.headshape   = a filename containing headshape, a structure containing a
%                     single triangulated boundary, or a Nx3 matrix with surface
%                     points
% If no headshape is specified, the dipole layer will be based on the inner compartment
% of the volume conduction model.
%
% Optionally, you can modify the leadfields by reducing the rank, i.e. remove the weakest orientation
%   cfg.reducerank    = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.backproject   = 'yes' or 'no',  determines when reducerank is applied whether the 
%                       lower rank leadfield is projected back onto the original linear 
%                       subspace, or not (default = 'yes')
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The following cfg fields are optional:
%   cfg.feedback
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_COMBINEPLANAR, FT_PREPARE_NEIGHBOURS

% Copyright (C) 2004-2019, Robert Oostenveld
% Copyright (C) 2020-      Robert Oostenveld and Jan-Mathijs Schoffelen
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
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% store the original input representation of the data, this is used later on to convert it back
isfreq = ft_datatype(data, 'freq');
istlck = ft_datatype(data, 'timelock');  % this will be temporary converted into raw

% check if the input data is valid for this function, this converts the data if needed
data = ft_checkdata(data, 'datatype', {'raw' 'freq'}, 'feedback', 'yes', 'hassampleinfo', 'yes', 'ismeg', 'yes', 'senstype', {'ctf151', 'ctf275', 'bti148', 'bti248', 'itab153', 'yokogawa160', 'yokogawa64'});

if isfreq
  if ~isfield(data, 'fourierspctrm'), ft_error('freq data should contain Fourier spectra'); end
end

cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'grid',    'sourcemodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'pruneratio', 'tolerance'});

% set the default configuration
cfg.channel      = ft_getopt(cfg, 'channel',      'MEG');
cfg.tolerance    = ft_getopt(cfg, 'tolerance',    1e-3);
cfg.trials       = ft_getopt(cfg, 'trials',       'all', 1);
cfg.planarmethod = ft_getopt(cfg, 'planarmethod', 'sincos');
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval',  {'headshape', 'headmodel', []});

if ~strcmp(cfg.planarmethod, 'sourceproject')
  % this is limited to reading neighbours from disk and/or selecting channels
  % the user should call FT_PREPARE_NEIGHBOURS directly for the actual construction
  tmpcfg = keepfields(cfg, {'neighbours', 'channel', 'showcallinfo'});
  cfg.neighbours = ft_prepare_neighbours(tmpcfg);
end

if isfield(cfg, 'headshape') && isa(cfg.headshape, 'config')
  % convert the nested config-object back into a normal structure
  cfg.headshape = struct(cfg.headshape);
end

if isfield(cfg, 'neighbours') && isa(cfg.neighbours, 'config')
  % convert the nested config-object back into a normal structure
  cfg.neighbours = struct(cfg.neighbours);
end

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'}); % this is moved to cfg.sourcemodel.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.sourcemodel.unit by the subsequent createsubcfg

% put the low-level options pertaining to the sourcemodel in their own field
cfg = ft_checkconfig(cfg, 'createsubcfg', {'sourcemodel'});
% move some fields from cfg.sourcemodel back to the top-level configuration
cfg = ft_checkconfig(cfg, 'createtopcfg', {'sourcemodel'});

% select trials of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo'}); % don't keep tolerance, it is used differently here
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

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
    ft_error('the method ''sourceproject'' is not supported for frequency data as input');
  end
  
  % PREPARE_HEADMODEL will match the data labels, the gradiometer labels and the
  % volume model labels (in case of a localspheres model) and result in a gradiometer
  % definition that only contains the gradiometers that are present in the data.
  [headmodel, axial.grad, cfg] = prepare_headmodel(cfg, data);
  
  % construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
  leadfieldopt = {};
  leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(cfg, 'reducerank'));
  leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(cfg, 'backproject'));
  leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(cfg, 'normalize'));
  leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(cfg, 'normalizeparam'));
  leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(cfg, 'weight'));

  % copy all options that are potentially used in FT_PREPARE_SOURCEMODEL
  tmpcfg           = keepfields(cfg, {'sourcemodel', 'mri', 'headshape', 'symmetry', 'smooth', 'threshold', 'spheremesh', 'inwardshift', 'xgrid' 'ygrid', 'zgrid', 'resolution', 'tight', 'warpmni', 'template', 'showcallinfo'});
  tmpcfg.headmodel = headmodel;
  tmpcfg.grad      = axial.grad;
  % determine the dipole layer that represents the surface of the brain
  sourcemodel = ft_prepare_sourcemodel(tmpcfg);
  
  % compute the forward model for the axial gradiometers
  ft_info('computing forward model for %d dipoles\n', size(sourcemodel.pos,1));
  lfold = ft_compute_leadfield(sourcemodel.pos, axial.grad, headmodel, leadfieldopt{:});
  
  % construct the planar gradient definition and compute its forward model
  % this will not work for a localspheres model, compute_leadfield will catch
  % the error
  planar.grad = constructplanargrad([], axial.grad);
  lfnew = ft_compute_leadfield(sourcemodel.pos, planar.grad, headmodel, leadfieldopt{:});
  
  % compute the interpolation matrix
  transform = lfnew * ft_inv(lfold, 'method', 'tsvd', 'tolerance', cfg.tolerance);
  
  planarmontage = [];
  planarmontage.tra      = transform;
  planarmontage.labelold = axial.grad.label;
  planarmontage.labelnew = planar.grad.label;
  
  % apply the linear transformation to the data
  interp  = ft_apply_montage(data, planarmontage, 'keepunused', 'yes');
  
  % also apply the linear transformation to the gradiometer definition
  interp.grad = ft_apply_montage(data.grad, planarmontage, 'balancename', 'planar', 'keepunused', 'yes');
  
  % ensure there is a type string describing the gradiometer definition
  if ~isfield(interp.grad, 'type')
    interp.grad.type = [ft_senstype(data.grad) '_planar'];
  else
    interp.grad.type = [interp.grad.type '_planar'];
  end
  
  %   % interpolate the data towards the planar gradiometers
  %   for i=1:Ntrials
  %     ft_info('interpolating trial %d to planar gradiometer\n', i);
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
  
  sens = ft_determine_units(data.grad);
  chanposnans = any(isnan(sens.chanpos(:))) || any(isnan(sens.chanori(:)));
  if chanposnans
    if isfield(sens, 'chanposold')
      % temporarily replace chanpos and chanorig with the original values
      sens.chanpos = sens.chanposold;
      sens.chanori = sens.chanoriold;
      sens.label   = sens.labelold;
      sens = rmfield(sens, {'chanposold', 'chanoriold', 'labelold'});
    else
      ft_error('The channel positions (and/or orientations) contain NaNs; this prohibits correct behavior of the function. Please replace the input channel definition with one that contains valid channel positions');
    end
  end
  cfg.channel = ft_channelselection(cfg.channel, sens.label);
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  
  % ensure channel order according to cfg.channel (there might be one check
  % too much in here somewhere or in the subfunctions, but I don't care.
  % Better one too much than one too little - JMH @ 09/19/12
  cfg = struct(cfg);
  [neighbsel] = match_str({cfg.neighbours.label}, cfg.channel);
  cfg.neighbours = cfg.neighbours(neighbsel);
  cfg.neighbsel = channelconnectivity(cfg);
  
  assert(any(cfg.neighbsel(:)), 'no neighbours found')
  ft_info('average number of neighbours is %.2f\n', mean(sum(cfg.neighbsel)));
  
  Ngrad = length(sens.label);
  distance = zeros(Ngrad,Ngrad);
  
  for i=1:size(cfg.neighbsel,1)
    j=find(cfg.neighbsel(i, :));
    d = sqrt(sum((sens.chanpos(j,:) - repmat(sens.chanpos(i, :), numel(j), 1)).^2, 2));
    distance(i,j) = d;
    distance(j,i) = d;
  end
  
  ft_info('minimum distance between neighbours is %6.2f %s\n', min(distance(distance~=0)), sens.unit);
  ft_info('maximum distance between gradiometers is %6.2f %s\n', max(distance(distance~=0)), sens.unit);
  
  % The following does not work when running in deployed mode because the
  % private functions that compute the planar montage are not recognized as
  % such and won't be compiled, unless explicitly specified.
  
  % % generically call megplanar_orig megplanar_sincos or megplanar_fitplane
  %fun = ['megplanar_'  cfg.planarmethod];
  %if ~exist(fun, 'file')
  %  ft_error('unknown method for computation of planar gradient');
  %end
  %planarmontage = eval([fun '(cfg, data.grad)']);
  
  switch cfg.planarmethod
    case 'sincos'
      planarmontage = megplanar_sincos(cfg, sens);
    case 'orig'
      % method specific info that is needed
      cfg.distance  = distance;
      
      planarmontage = megplanar_orig(cfg, sens);
    case 'fitplane'
      planarmontage = megplanar_fitplane(cfg, sens);
    otherwise
      fun = ['megplanar_' cfg.planarmethod];
      if ~exist(fun, 'file')
        ft_error('unknown method for computation of planar gradient');
      end
      planarmontage = eval([fun '(cfg, data.grad)']);
  end
  
  % apply the linear transformation to the data
  interp = ft_apply_montage(data, planarmontage, 'keepunused', 'yes', 'feedback', cfg.feedback);
  
  % also apply the linear transformation to the gradiometer definition
  interp.grad = ft_apply_montage(sens, planarmontage, 'balancename', 'planar', 'keepunused', 'yes');
  
  % ensure there is a type string describing the gradiometer definition
  if ~isfield(interp.grad, 'type')
    % put the original gradiometer type in (will get _planar appended)
    interp.grad.type = ft_senstype(sens);
  end
  interp.grad.type = [interp.grad.type '_planar'];
  
  % add the chanpos info back into the gradiometer description
  tmplabel = interp.grad.label;
  for k = 1:numel(tmplabel)
    if ~isempty(strfind(tmplabel{k}, '_dV')) || ~isempty(strfind(tmplabel{k}, '_dH'))
      tmplabel{k} = tmplabel{k}(1:end-3);
    end
  end
  [ix,iy] = match_str(tmplabel, sens.label);
  interp.grad.chanpos(ix,:) = sens.chanpos(iy,:);
  
  % if the original chanpos contained nans, make sure to put nans in the
  % updated one as well, and move the updated chanpos values to chanposold
  if chanposnans
    interp.grad.chanposold = sens.chanpos;
    interp.grad.chanoriold = sens.chanori;
    interp.grad.labelold   = sens.label;
    interp.grad.chanpos    = nan(size(interp.grad.chanpos));
    interp.grad.chanori    = nan(size(interp.grad.chanori));
  end
end

if istlck
  % convert the raw structure back into a timelock structure
  interp = ft_checkdata(interp, 'datatype', 'timelock');
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
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = interp;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
