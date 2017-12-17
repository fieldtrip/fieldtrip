function [data] = ft_channelrepair(cfg, data)

% FT_CHANNELREPAIR repairs bad or missing channels in the data by replacing them with the
% plain average of of all neighbours, by a weighted average of all neighbours, by an
% interpolation based on a surface Laplacian, or by spherical spline interpolating (see
% Perrin et al., 1989).
%
% Use as
%   [interp] = ft_channelrepair(cfg, data)
%
% The configuration must contain
%   cfg.method         = 'weighted', 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
%   cfg.badchannel     = cell-array, see FT_CHANNELSELECTION for details
%   cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
%   cfg.neighbours     = neighbourhood structure, see also FT_PREPARE_NEIGHBOURS
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.lambda         = regularisation parameter (default = 1e-5, not for method 'distance')
%   cfg.order          = order of the polynomial interpolation (default = 4, not for method 'distance')
%
% The weighted neighbour approach cannot be used reliably to repair multiple bad channels
% that lie next to each other.
%
% If you want to reconstruct channels that are absent in your data, those channels may
% also be missing from the sensor definition (grad, elec or opto) and determining the
% neighbours is non-trivial. In that case you must use a complete sensor definition from
% another dataset or from a template.
%
% The EEG, MEG or NIRS sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
%   cfg.opto          = structure with optode definition, see FT_DATATYPE_SENS
%   cfg.optofile      = name of file containing the optode definition, see FT_READ_SENS
%
% This function only interpolates data over space, not over time. If you want to
% interpolate using temporal information, e.g. using a segment of data before and
% after the nan-marked artifact, you should use FT_INTERPOLATENAN.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_MEGREALIGN, FT_MEGPLANAR, FT_PREPARE_NEIGHBOURS, FT_INTERPOLATENAN

% Copyright (C) 2004-2009, Robert Oostenveld
% Copyright (C) 2012-2013, J?rn M. Horschig, Jason Farquhar

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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'nearest', 'weighted'});

% set the default configuration
cfg.badchannel     = ft_getopt(cfg, 'badchannel',     {});
cfg.missingchannel = ft_getopt(cfg, 'missingchannel', {});
cfg.trials         = ft_getopt(cfg, 'trials',         'all', 1);
cfg.method         = ft_getopt(cfg, 'method',         'weighted');
cfg.lambda         = ft_getopt(cfg, 'lambda',         []); % subfunction will handle this
cfg.order          = ft_getopt(cfg, 'order',          []); % subfunction will handle this

% check if the input cfg is valid for this function
if strcmp(cfg.method, 'weighted')
  cfg = ft_checkconfig(cfg, 'required', {'neighbours'});
end

% store the original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% select trials of interest
tmpcfg = [];
tmpcfg.trials = cfg.trials;
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

if strcmp(cfg.method, 'nan')
  % this does not require the spatial information of the channels
  sens = [];
  sens.chanpos = zeros(0, 3);
  sens.label = cell(0, 1);
else
  % this requires the spatial information of the channels
  try
    % prefer sens from cfg over sens from data
    sens = ft_fetch_sens(cfg);
  catch
    sens = ft_fetch_sens(cfg, data);
  end
end

% determine the type of data
iseeg  = ft_senstype(sens, 'eeg');
ismeg  = ft_senstype(sens, 'meg');
isnirs = ft_senstype(sens, 'opto');

% check if any of the channel positions contains NaNs; this happens when
% component data are backprojected to the sensor level
if any(isnan(sens.chanpos(:)))
  ft_error('The channel positions contain NaNs; this prohibits correct behavior of the function. Please replace the input channel definition with one that contains valid channel positions');
end

if ismeg && ~any(strcmp(ft_senstype(sens), {'ctf151', 'ctf275', 'bti148', 'bti248', 'babysquid74'}))
  % MEG systems with only magnetometers or axial gradiometers are easy, planar systems are not
  ft_warning('be careful when using "%s" - mixing of sensor types (e.g. magnetometers and gradiometers) can lead to wrong data. Check your neighbour-structure thoroughly', ft_senstype(sens));
end

% get selection of channels that are missing and/or bad
cfg.missingchannel = cat(1, cfg.missingchannel(:), cfg.badchannel);
cfg.missingchannel = setdiff(cfg.missingchannel, data.label);
cfg.badchannel = ft_channelselection(cfg.badchannel, data.label);
fprintf('There are %d bad channels\n', length(cfg.badchannel));
fprintf('There are %d missing channels\n', length(cfg.missingchannel));

% warn if weighted neighbour approach (see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=634)
if ~isempty(cfg.missingchannel) && strcmp(cfg.method, 'weighted')
  ft_warning('Reconstructing missing channels using the weighted neighbour approach is not recommended!');
end

% store the realigned data in a new structure
interp = [];
interp.label   = data.label;
interp.time    = data.time;

% first repair badchannels
if strcmp(cfg.method, 'weighted') || strcmp(cfg.method, 'average')
  
  if ~isempty(cfg.badchannel)
    [goodchanlabels, goodchanindcs] = setdiff(data.label, cfg.badchannel);
    goodchanindcs = sort(goodchanindcs); % undo automatical sorting by setdiff
    connectivityMatrix = channelconnectivity(cfg, data);
    connectivityMatrix = connectivityMatrix(:, goodchanindcs); % all chans x good chans
    
    ntrials = length(data.trial);
    nchans  = length(data.label);
    
    repair  = eye(nchans, nchans);
    badindx = match_str(data.label, cfg.badchannel);
    
    for k=badindx'
      fprintf('repairing channel %s\n', data.label{k});
      repair(k, k) = 0;
      l = goodchanindcs(connectivityMatrix(k, :));
      % get bad channels out
      [a, b] = setdiff(data.label(l), data.label(badindx));
      b = sort(b); % undo automatical sorting by setdiff
      l(~ismember(find(l), b)) = [];
      % get corresponding ids for sens structure
      [a, b] = match_str(sens.label, data.label(l));
      goodsensindx = a(b);
      [a, b] = match_str(sens.label, data.label(k));
      badsensindx = a(b);
      fprintf('\tusing neighbour %s\n', sens.label{goodsensindx});
      if strcmp(cfg.method, 'weighted')
        distance = sqrt(sum((sens.chanpos(goodsensindx, :) - repmat(sens.chanpos(badsensindx, :), numel(goodsensindx), 1)).^2, 2));
      elseif strcmp(cfg.method, 'average')
        distance = 1;
      end
      repair(k, l) = (1./distance);
      repair(k, l) = repair(k, l) ./ sum(repair(k, l));
    end
    
    % use sparse matrix to speed up computations
    repair = sparse(repair);
    
    % compute the repaired data for each trial
    fprintf('\n');
    fprintf('repairing bad channels for %i trials %d', ntrials);
    for i=1:ntrials
      fprintf('.');
      interp.trial{i} = repair * data.trial{i};
    end
    fprintf('\n');
  else
    fprintf('no bad channels to repair\n');
    interp.trial = data.trial;
  end
  
  if ~isempty(cfg.missingchannel)
    fprintf('Interpolated missing channels will be concatenated.\n');
    
    nchans  = length(interp.label);
    ntrials = length(interp.trial);
    
    % interpolation missing channels
    goodchanindcs = 1:numel(data.label);
    for chan=1:numel(cfg.missingchannel)
      interp.label{end+1} = cfg.missingchannel{chan};
      % creating dummy trial data
      for i=1:ntrials
        interp.trial{i}(end+1, :) = 0;
      end
    end
    connectivityMatrix = channelconnectivity(cfg, interp);
    connectivityMatrix = connectivityMatrix(:, goodchanindcs); % all chans x good chans
    
    repair  = eye(nchans, nchans);
    missingindx = match_str(interp.label, cfg.missingchannel);
    unable = [];
    for k=missingindx'
      fprintf('trying to reconstruct missing channel %s\n', interp.label{k});
      repair(k, k) = 0;
      l = goodchanindcs(connectivityMatrix(k, :));
      % get bad channels out
      [a, b] = setdiff(data.label(l), interp.label(missingindx));
      b = sort(b); % undo automatical sorting by setdiff
      l(~ismember(find(l), b)) = [];
      % get corresponding ids for sens structure
      [a, b] = match_str(sens.label, interp.label(l));
      goodsensindx = a(b);
      if isempty(goodsensindx)
        fprintf('\tcannot reconstruct channel - no neighbours in the original data or in the sensor position\n');
        unable = [unable k];
      else
        [a, b] = match_str(sens.label, interp.label(k));
        badsensindx = a(b);
        fprintf('\tusing neighbour %s\n', sens.label{goodsensindx});
        if strcmp(cfg.method, 'weighted')
          distance = sqrt(sum((sens.chanpos(goodsensindx, :) - repmat(sens.chanpos(badsensindx, :), numel(goodsensindx), 1)).^2, 2));
        elseif strcmp(cfg.method, 'average')
          distance = 1;
        end
        repair(k, l) = (1./distance);
        repair(k, l) = repair(k, l) ./ sum(repair(k, l));
      end
    end
    
    % use sparse matrix to speed up computations
    repair = sparse(repair);
    
    fprintf('\n');
    % compute the missing data for each trial and remove those could not be
    % reconstructed
    fprintf('\n');
    fprintf('interpolating missing channel for %i trials %d', ntrials);
    for i=1:ntrials
      fprintf('.');
      interp.trial{i} = repair * interp.trial{i};
      interp.trial{i}(unable, :) = [];
    end
    
    interp.label(unable) = [];
    fprintf('\n');
  end
  
elseif strcmp(cfg.method, 'spline') || strcmp(cfg.method, 'slap')
  if ~isempty(cfg.badchannel) || ~isempty(cfg.missingchannel)
    fprintf('Spherical spline and surface Laplacian interpolation will treat bad and missing channels the same. Missing channels will be concatenated at the end of your data structure.\n');
  end
  % subselect only those sensors that are in the data or in badchannel or missingchannel
  badchannels   = union(cfg.badchannel, cfg.missingchannel);
  sensidx       = ismember(sens.label, union(data.label, badchannels));
  label    = sens.label(sensidx);
  chanpos  = sens.chanpos(sensidx, :);
  try, chanori   = sens.chanori(sensidx, :); end
  try, chantype  = sens.chantype(sensidx, :); end
  try, chanunit  = sens.chanunit(sensidx, :); end
  
  fprintf('Checking spherical fit... ');
  [c, r] = fitsphere(chanpos);
  d = chanpos - repmat(c, numel(find(sensidx)), 1);
  d = sqrt(sum(d.^2, 2));
  d = mean(abs(d) / r);
  if abs(d-1) > 0.1
    ft_warning('bad spherical fit (residual: %.2f%%). The interpolation will be inaccurate.', 100*(d-1));
  elseif abs(d-1) < 0.01
    fprintf('perfect spherical fit (residual: %.1f%%)\n', 100*(d-1));
  else
    fprintf('good spherical fit (residual: %.1f%%)\n', 100*(d-1));
  end
  
  if strcmp(cfg.method, 'slap')
    ft_warning('''slap'' method is not fully supported - be careful in interpreting your results');
  end
  % move missing channels to the end
  missidx = find(ismember(label, cfg.missingchannel));
  label(end+1:end+numel(missidx))      = label(missidx);
  label(missidx)                       = [];
  chanpos(end+1:end+numel(missidx), :) = chanpos(missidx, :);
  chanpos(missidx, :)                  = [];
  
  % select good channels only for interpolation
  [goodchanlabels, goodchanindcs] = setdiff(label, badchannels);
  allchans = false;
  if isempty(goodchanindcs)
    goodchanindcs = 1:numel(label);
    allchans = true;
    ft_warning('No good channels found - interpolating based on all channels');
  end
  % undo automatical sorting by setdiff
  goodchanindcs      = sort(goodchanindcs);
  % only take good channels that are in data (and remember how they are sorted)
  [dataidx, sensidx] = match_str(data.label, label(goodchanindcs));
  
  % interpolate
  fprintf('computing weight matrix...');
  repair = sphericalSplineInterpolate(chanpos(goodchanindcs(sensidx), :)', chanpos', cfg.lambda, cfg.order, cfg.method);
  fprintf(' done!\n');
  
  if ~allchans
    % only use the rows corresponding to the channels that actually need interpolation
    repair(goodchanindcs(sensidx), :) = 0;
    for k = 1:numel(sensidx)
      i = strcmp(label(goodchanindcs(sensidx(k))), label(goodchanindcs(sensidx)));
      repair(goodchanindcs(sensidx(k)), i) = 1;
    end
  end % else all rows need to be interpolated
  
  % compute the missing data for each trial and remove those could not be reconstructed
  ntrials = length(data.trial);
  fprintf('\n');
  fprintf('interpolating channels for %i trials %d', ntrials);
  for i=1:ntrials
    fprintf('.');
    interp.trial{i} = repair * data.trial{i}(dataidx, :);
  end
  fprintf('\n');
  % update channels labels due to reordering by
  interp.label = label;
  
elseif strcmp(cfg.method, 'nan')
  % copy the original data over to the output data structure
  interp.trial = data.trial;
  if ~isempty(cfg.badchannel)
    for k=1:length(cfg.badchannel)
      fprintf('inserting nan into bad channel %s\n', cfg.badchannel{k});
    end
    badindx = match_str(interp.label, cfg.badchannel);
    for i=1:numel(interp.trial)
      interp.trial{i}(badindx, :) = nan;
    end
  end % bad channels
  if ~isempty(cfg.missingchannel)
    for k=1:length(cfg.missingchannel)
      fprintf('inserting nan into missing channel %s\n', cfg.missingchannel{k});
    end
    interp.label = cat(1, interp.label(:), cfg.missingchannel(:));
    for i=1:numel(interp.trial)
      interp.trial{i} = cat(1, interp.trial{i}, nan(length(cfg.missingchannel), length(interp.time{i})));
    end
  end % missing channels
  
else
  ft_error('unknown method "%s" for interpolation', cfg.method);
end

% copy the additional fields over to the newly interpolated data
datafields   = fieldnames(data);
interpfields = fieldnames(interp);
exfields     = setdiff(datafields, interpfields);
for f = 1:length(exfields)
  interp.(exfields{f}) = data.(exfields{f});
end

% re-insert the sensor array
if iseeg
  interp.elec  = sens;
elseif ismeg
  interp.grad  = sens;
elseif isnirs
  interp.opto  = sens;
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    interp = ft_checkdata(interp, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
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
