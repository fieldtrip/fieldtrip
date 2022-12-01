function [data] = ft_channelrepair(cfg, data)

% FT_CHANNELREPAIR repairs bad or missing channels in the data by replacing
% them with the plain average of of all neighbours, by a weighted average
% of all neighbours, by an interpolation based on a surface Laplacian, or
% by spherical spline interpolating (see Perrin et al., 1989).
%
% Use as
%   [interp] = ft_channelrepair(cfg, data)
% where the input data corresponds to the output from FT_PREPROCESSING.
%
% The configuration should contain
%   cfg.method         = 'weighted', 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
%   cfg.badchannel     = cell-array, see FT_CHANNELSELECTION for details
%   cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
%   cfg.neighbours     = neighbourhood structure, see FT_PREPARE_NEIGHBOURS for details
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.lambda         = regularisation parameter (default = 1e-5, for method 'spline' and 'slap')
%   cfg.order          = order of the polynomial interpolation (default = 4 for methods 'spline' and 'slap')
%   cfg.senstype       = string, which type of data to repair. Can be 'meg', 'eeg' or 'nirs' (default is automatic)
%
% The weighted neighbour approach cannot be used reliably to repair multiple
% bad channels that lie next to each other.
%
% If you want to reconstruct channels that are absent in your data, those
% channels may also be missing from the sensor definition (grad, elec or opto)
% and determining the neighbours is non-trivial. In that case you must use
% a complete sensor definition from another dataset or from a template.
%
% The EEG, MEG or NIRS sensor positions can be present as a field in the
% data (data.grad/data.elec/data.opto, depending on the type of data),
% or can be specified as cfg option. Either one is required for the following
% methods: 'weighted', 'spline', and 'slap'. Depending on the type of
% data this should be one of the following
%   cfg.elec = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.grad = structure with gradiometer definition or filename, see FT_READ_SENS
%   cfg.opto = structure with optode definition, see FT_READ_SENS
%
% This function will only repair one type of channels (MEG, EEG or NIRS) at
% a time. If you want to repair multiple types of channels, you should call
% it multiple times and use FT_SELECTDATA and FT_APPENDDATA.
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
% Copyright (C) 2012-2013, Jörn M. Horschig, Jason Farquhar
% Copyright (C) 2021, Jan-Mathijs Schoffelen

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

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'nearest', 'weighted'});
cfg = ft_checkconfig(cfg, 'renamed',    {'elecfile', 'elec'});
cfg = ft_checkconfig(cfg, 'renamed',    {'gradfile', 'grad'});
cfg = ft_checkconfig(cfg, 'renamed',    {'optofile', 'opto'});

% set the default configuration
cfg.badchannel     = ft_getopt(cfg, 'badchannel',     {});
cfg.missingchannel = ft_getopt(cfg, 'missingchannel', {});
cfg.senstype       = ft_getopt(cfg, 'senstype',       []); % default is handled below
cfg.trials         = ft_getopt(cfg, 'trials',         'all', 1);
cfg.method         = ft_getopt(cfg, 'method',         'weighted');
cfg.lambda         = ft_getopt(cfg, 'lambda',         []); % subfunction will handle this
cfg.order          = ft_getopt(cfg, 'order',          []); % subfunction will handle this

% check if the input cfg is valid for this function
if ismember(cfg.method, {'weighted', 'average'})
  cfg = ft_checkconfig(cfg, 'required', {'neighbours'});
end

cfg = ft_checkconfig(cfg, 'forbidden', 'layout');

% store the original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% select trials of interest
tmpcfg = keepfields(cfg, {'trials', 'tolerance', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

if isempty(cfg.badchannel)
  % check if the first sample of the first trial contains NaNs; if so treat it as a bad channel
  cfg.badchannel = ft_channelselection(find(isnan(data.trial{1}(:,1))), data.label);
  if ~isempty(cfg.badchannel)
    ft_info('detected channel %s as bad\n', cfg.badchannel{:});
  end
end

needsens = true;
if ismember(cfg.method, {'nan', 'average'})
  % a sensor description is not needed
  needsens = false;
end

if needsens
  
  % sometimes the data can have both gradiometers, electrodes and/or optodes
  % but this function can only deal with one type of data at a time
  if isempty(cfg.senstype)
    % look at the bad channels, then at the data
    if all(ft_chantype(cat(1, cfg.missingchannel(:), cfg.badchannel(:)), 'meg')) || ft_senstype(data, 'meg')
      ft_notice('assuming that the data is MEG, use cfg.senstype to overrule this');
      cfg.senstype = 'meg';
    elseif all(ft_chantype(cat(1, cfg.missingchannel(:), cfg.badchannel(:)), 'eeg')) || ft_senstype(data, 'eeg')
      ft_notice('assuming that the data is EEG, use cfg.senstype to overrule this');
      cfg.senstype = 'eeg';
    elseif all(ft_chantype(cat(1, cfg.missingchannel(:), cfg.badchannel(:)), 'nirs')) || ft_senstype(data, 'nirs')
      ft_notice('assuming that the data is NIRS, use cfg.senstype to overrule this');
      cfg.senstype = 'nirs';
    else
      % let FT_FETCH_SENS decide which sens to return
    end
  end
  
  % the 3D spatial information of the channels is needed
  % this will prefer sens from cfg over sens from data
  sens = ft_fetch_sens(cfg, data);
  
  % check if any of the channel positions contains NaNs; this happens when
  % component data are backprojected to the sensor level
  if any(isnan(sens.chanpos(:)))
    ft_error('The channel positions contain NaNs; this prohibits correct behavior of the function. Please replace the input channel definition with one that contains valid channel positions');
  end
  
  % determine the global type of the data
  iseeg  = ft_senstype(sens, 'eeg');
  ismeg  = ft_senstype(sens, 'meg');
  isnirs = ft_senstype(sens, 'opto');
  
  % determine the detailled type, e.g. ctf151 or neuromag122
  sensortype = ft_senstype(sens);
  
else
  % determine the global type of the data
  iseeg  = ft_senstype(data, 'eeg');
  ismeg  = ft_senstype(data, 'meg');
  isnirs = ft_senstype(data, 'opto');
  % determine the detailled type, e.g. ctf151 or neuromag122
  sensortype = ft_senstype(data);
end

if ismeg && ~any(strcmp(sensortype, {'ctf151', 'ctf275', 'bti148', 'bti248', 'babysquid74'}))
  % MEG systems with only magnetometers or axial gradiometers are easy, planar systems are not
  ft_warning('be careful when using "%s" - mixing of sensor types (e.g. magnetometers and gradiometers) can lead to wrong data. Check your neighbour-structure thoroughly', ft_senstype(sens));
end

if ~isempty(cfg.missingchannel) && strcmp(cfg.method, 'weighted')
  ft_warning('Reconstructing missing channels using the weighted neighbour approach is not recommended!');
end

% get selection of channels that are missing and/or bad
cfg.missingchannel = cat(1, cfg.missingchannel(:), cfg.badchannel(:));
cfg.missingchannel = setdiff(cfg.missingchannel, data.label, 'stable');
cfg.badchannel     = ft_channelselection(cfg.badchannel, data.label);
fprintf('There are %d bad channels\n',     length(cfg.badchannel));
fprintf('There are %d missing channels\n', length(cfg.missingchannel));

% warn if weighted neighbour approach (see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=634)
if ~isempty(cfg.missingchannel) && strcmp(cfg.method, 'weighted')
  ft_warning('Reconstructing missing channels using the weighted neighbour approach is not recommended!');
end

% store the realigned data in a new structure
interp         = [];
interp.label   = data.label;
interp.time    = data.time;

switch cfg.method
  case {'weighted' 'average'}
    
    % first repair badchannels
    if ~isempty(cfg.badchannel)
      [dum, goodindx] = setdiff(data.label, cfg.badchannel, 'stable');
      adj = channelconnectivity(cfg, data);
      adj = adj(:, goodindx); % all chans x good chans
      
      ntrials = length(data.trial);
      nchans  = length(data.label);
      
      repair  = eye(nchans, nchans);
      badindx = match_str(data.label, cfg.badchannel);
      unable  = zeros(1,0);
      
      for k = badindx(:)'
        fprintf('repairing channel %s\n', data.label{k});
        repair(k, k) = 0;
        
        list = goodindx(adj(k, :));
        % get bad channels out
        [a, b] = setdiff(data.label(list), data.label(badindx));
        b      = sort(b); % undo automatical sorting by setdiff
        list(~ismember(find(list), b)) = [];
        
        distance = 1; % needed here, in case isempty(goodsensindx)
        if strcmp(cfg.method, 'weighted')
          % get corresponding ids for sens structure
          [a, b] = match_str(sens.label, data.label(list));
          goodsensindx = a(b);
          if isempty(goodsensindx)
            fprintf('\tcannot reconstruct channel - no neighbours in the original data or in the sensor positions\n');
            unable = cat(2, unable, k);
          else
            [a, b] = match_str(sens.label, data.label(k));
            badsensindx = a(b);
            fprintf('\tusing neighbour %s\n', sens.label{goodsensindx});
            distance = sqrt(sum((sens.chanpos(goodsensindx, :) - repmat(sens.chanpos(badsensindx, :), numel(goodsensindx), 1)).^2, 2));
          end
        end
        repair(k, list) = (1./distance);
        repair(k, list) = repair(k, list) ./ sum(repair(k, list));
      end
      
      % use sparse matrix to speed up computations
      repair = sparse(repair);
      
      % compute the missing data for each trial and set the ones that could not be reconstructed to nan
      fprintf('\n');
      fprintf('interpolating bad channels for %i trials %d', ntrials);
      for k = 1:ntrials
        fprintf('.');
        interp.trial{k} = repair * data.trial{k};
        interp.trial{k}(unable, :) = nan;
      end
      fprintf('\n');
    else
      fprintf('no bad channels to repair\n');
      interp.trial = data.trial;
    end
    
    % deal with missingchannels next
    if ~isempty(cfg.missingchannel)
      fprintf('Interpolated missing channels will be concatenated.\n');
      
      nchans  = length(interp.label);
      ntrials = length(interp.trial);
      
      % interpolation missing channels
      goodindx = 1:numel(data.label);
      for chan = 1:numel(cfg.missingchannel)
        interp.label{end+1} = cfg.missingchannel{chan};
        % creating dummy trial data
        for k = 1:ntrials
          interp.trial{k}(end+1, :) = 0;
        end
      end
      adj = channelconnectivity(cfg, interp);
      adj = adj(:, goodindx); % all chans x good chans
      
      repair      = eye(nchans, nchans);
      missingindx = match_str(interp.label, cfg.missingchannel);
      unable      = zeros(1,0);
      
      for k = missingindx(:)'
        fprintf('trying to reconstruct missing channel %s\n', interp.label{k});
        repair(k, k) = 0;
        list = goodindx(adj(k, :));
        % get bad channels out
        [a, b] = setdiff(data.label(list), interp.label(missingindx));
        b      = sort(b); % undo automatical sorting by setdiff
        list(~ismember(find(list), b)) = [];
        
        distance = 1; % needed here, in case isempty(goodsensindx)
        if strcmp(cfg.method, 'weighted')
          % get corresponding ids for sens structure
          [a, b] = match_str(sens.label, interp.label(list));
          goodsensindx = a(b);
          if isempty(goodsensindx)
            fprintf('\tcannot reconstruct channel - no neighbours in the original data or in the sensor positions\n');
            unable = cat(2, unable, k);
          else
            [a, b] = match_str(sens.label, interp.label(k));
            badsensindx = a(b);
            fprintf('\tusing neighbour %s\n', sens.label{goodsensindx});
            distance = sqrt(sum((sens.chanpos(goodsensindx, :) - repmat(sens.chanpos(badsensindx, :), numel(goodsensindx), 1)).^2, 2));
          end
        end
        repair(k, list) = (1./distance);
        repair(k, list) = repair(k, list) ./ sum(repair(k, list));
      end
      
      % use sparse matrix to speed up computations
      repair = sparse(repair);
      
      fprintf('\n');
      % compute the missing data for each trial and set the ones that could not be reconstructed to nan
      fprintf('\n');
      fprintf('interpolating missing channel for %i trials %d', ntrials);
      for k = 1:ntrials
        fprintf('.');
        interp.trial{k} = repair * interp.trial{k};
        interp.trial{k}(unable, :) = nan;
      end
      fprintf('\n');
    end
    
  case {'spline' 'slap'}
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
    [goodchanlabels, goodindx] = setdiff(label, badchannels);
    allchans = false;
    if isempty(goodindx)
      goodindx = 1:numel(label);
      allchans = true;
      ft_warning('No good channels found - interpolating based on all channels');
    end
    % undo automatical sorting by setdiff
    goodindx      = sort(goodindx);
    % only take good channels that are in data (and remember how they are sorted)
    [dataidx, sensidx] = match_str(data.label, label(goodindx));
    
    % interpolate
    fprintf('computing weight matrix...');
    [pot,lap] = sphsplint(chanpos(goodindx(sensidx), :), chanpos, cfg.order, 500, cfg.lambda);
    fprintf(' done!\n');
    
    % Chooses the right output.
    if strcmp(cfg.method, 'spline'), repair = pot;
    else, repair = lap;
    end
    
    if ~allchans
      % only use the rows corresponding to the channels that actually need interpolation
      repair(goodindx(sensidx), :) = 0;
      for k = 1:numel(sensidx)
        i = strcmp(label(goodindx(sensidx(k))), label(goodindx(sensidx)));
        repair(goodindx(sensidx(k)), i) = 1;
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
    % update channels labels due to reordering
    interp.label = label;
    
  case 'nan'
    % copy the original data over to the output data structure
    interp.trial = data.trial;
    if ~isempty(cfg.badchannel)
      for k = 1:length(cfg.badchannel)
        fprintf('inserting nan into bad channel %s\n', cfg.badchannel{k});
      end
      badindx = match_str(interp.label, cfg.badchannel);
      for k = 1:numel(interp.trial)
        interp.trial{k}(badindx, :) = nan;
      end
    end % bad channels
    if ~isempty(cfg.missingchannel)
      for k = 1:length(cfg.missingchannel)
        fprintf('inserting nan into missing channel %s\n', cfg.missingchannel{k});
      end
      interp.label = cat(1, interp.label(:), cfg.missingchannel(:));
      for k = 1:numel(interp.trial)
        interp.trial{k} = cat(1, interp.trial{k}, nan(length(cfg.missingchannel), length(interp.time{k})));
      end
    end % missing channels
    
  otherwise
    ft_error('unknown method "%s" for interpolation', cfg.method);
end

% copy the additional fields over to the newly interpolated data
datafields   = fieldnames(data);
interpfields = fieldnames(interp);
exfields     = setdiff(datafields, interpfields);
for k = 1:length(exfields)
  interp.(exfields{k}) = data.(exfields{k});
end

if needsens
  % re-insert the sensor array
  if iseeg
    interp.elec  = sens;
  elseif ismeg
    interp.grad  = sens;
  elseif isnirs
    interp.opto  = sens;
  end
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
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = interp;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
