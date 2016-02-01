function [data] = ft_rejectcomponent(cfg, comp, data)

% FT_REJECTCOMPONENT backprojects an ICA (or similar) decomposition to the
% channel level after removing the independent components that contain
% the artifacts. This function does not automatically detect the artifact
% components, you will have to do that yourself.
%
% Use as
%    [data] = ft_rejectcomponent(cfg, comp)
% or as
%    [data] = ft_rejectcomponent(cfg, comp, data)
%
% where the input comp is the result of FT_COMPONENTANALYSIS. The output
% data will have the same format as the output of FT_PREPROCESSING.
%
% An optional input argument data can be provided. In that case
% componentanalysis will do a subspace projection of the input data
% onto the space which is spanned by the topographies in the unmixing
% matrix in comp, after removal of the artifact components.  Please use
% this option of including data as input, if you wish to use the output
% data.grad in further computation, for example for leadfield computation.
%
% The configuration structure can contain
%   cfg.component  = list of components to remove, e.g. [1 4 7] or see FT_CHANNELSELECTION
%   cfg.demean     = 'no' or 'yes', whether to demean the input data (default = 'yes')
%   cfg.updatesens = 'no' or 'yes' (default = 'yes')
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_COMPONENTANALYSIS, FT_PREPROCESSING

% Copyright (C) 2005-2014, Robert Oostenveld
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
ft_preamble init
ft_preamble debug
ft_preamble loadvar comp data
ft_preamble provenance comp data
ft_preamble trackconfig

% the abort variable is set to true or false in ft_preamble_init
if abort
  return
end

% set the defaults
cfg.component       = ft_getopt(cfg, 'component',  []);
cfg.demean          = ft_getopt(cfg, 'demean',    'yes');
cfg.feedback        = ft_getopt(cfg, 'feedback',  'text');
cfg.updatesens      = ft_getopt(cfg, 'updatesens',  'yes');

% the data can be passed as input arguments or can be read from disk
nargin = 1;
nargin = nargin + exist('comp', 'var');
nargin = nargin + exist('data', 'var');


if nargin==3
  % check if the input data is valid for this function
  istlck  = ft_datatype(data, 'timelock');  % this will be temporary converted into raw
  data    = ft_checkdata(data, 'datatype', 'raw');
  comp    = ft_checkdata(comp, 'datatype', 'comp');
  label   = data.label;
  nchans  = length(data.label);
  ncomps  = length(comp.label);
  hasdata = 1;
elseif nargin==2
  % check if the input data is valid for this function
  istlck  = ft_datatype(comp, 'timelock');  % this will be temporary converted into raw
  comp    = ft_checkdata(comp, 'datatype', 'raw+comp');
  label   = comp.topolabel;
  ncomps  = length(comp.label);
  hasdata = 0;
else
  error('incorrect number of input arguments');
end

% cfg.component can be indicated by number or by label
cfg.component = ft_channelselection(cfg.component, comp.label);
reject = match_str(comp.label, cfg.component);

if isempty(reject)
  warning('no components were selected for rejection');
end

if min(reject)<1
  error('you cannot remove components that are not present in the data');
end

if max(reject)>ncomps
  error('you cannot remove components that are not present in the data');
end

if nargin==3 && strcmp(cfg.demean, 'yes')
  % optionally perform baseline correction on each trial
  fprintf('baseline correcting data \n');
  for trial=1:numel(data.trial)
    data.trial{trial} = ft_preproc_baselinecorrect(data.trial{trial});
  end
end

% set the rejected component amplitudes to zero
fprintf('removing %d components\n', length(reject));
if ~hasdata,
  fprintf('keeping %d components\n',  ncomps-length(reject));
else
  fprintf('keeping %d components\n',  nchans-length(reject));
end

% create a projection matrix by subtracting the subspace spanned by the
% topographies of the to-be-removed components from identity
[seldat, selcomp] = match_str(label, comp.topolabel);

if length(seldat)~=length(label) && nargin==3,
  warning('the subspace projection is not guaranteed to be correct for non-orthogonal components');
end

if hasdata
  mixing   = comp.topo(selcomp,:);
  unmixing = comp.unmixing(:,selcomp);
  
  % I am not sure about this, but it gives comparable results to the ~hasdata case
  % when comp contains non-orthogonal (=ica) topographies, and contains a complete decomposition
  
  montage     = [];
  montage.tra = eye(length(selcomp)) - mixing(:, reject)*unmixing(reject, :);
  % we are going from data to components, and back again
  montage.labelorg = comp.topolabel(selcomp);
  montage.labelnew = comp.topolabel(selcomp);
  
  keepunused = 'yes'; % keep the original data which are not present in the mixing provided
  
else
  mixing = comp.topo(selcomp, :);
  mixing(:, reject) = 0;
  
  montage     = [];
  montage.tra = mixing;
  % we are going from components to data
  montage.labelorg = comp.label;
  montage.labelnew = comp.topolabel(selcomp);
  
  keepunused = 'no'; % don't need to keep the original rejected components
  
  % create data structure
  data         = [];
  data.trial   = comp.trial;
  data.time    = comp.time;
  data.label   = comp.label;
  data.fsample = comp.fsample;
  if isfield(comp, 'grad'),       data.grad       = comp.grad;       end
  if isfield(comp, 'elec'),       data.elec       = comp.elec;       end
  if isfield(comp, 'trialinfo'),  data.trialinfo  = comp.trialinfo;  end
  if isfield(comp, 'sampleinfo'), data.sampleinfo = comp.sampleinfo; end
end % if hasdata

% apply the montage to the data and to the sensor description
data = ft_apply_montage(data, montage, 'keepunused', keepunused, 'feedback', cfg.feedback);

% apply the montage also to the elec/grad, if present
if isfield(data, 'grad')
  sensfield = 'grad';
  if strcmp(cfg.updatesens, 'yes')
    fprintf('applying the backprojection matrix to the gradiometer description\n');
  else
    fprintf('not applying the backprojection matrix to the gradiometer description\n');
  end
elseif isfield(data, 'elec') && isfield(data.elec, 'tra')
  sensfield = 'elec';
  if strcmp(cfg.updatesens, 'yes')
    fprintf('applying the backprojection matrix to the electrode description\n');
  else
    fprintf('not applying the backprojection matrix to the electrode description\n');
  end
else
  fprintf('not applying the backprojection matrix to the sensor description\n');
  sensfield = [];
end

if ~isempty(sensfield) && strcmp(cfg.updatesens, 'yes')
  % keepunused = 'yes' is required to get back e.g. reference or otherwise
  % unused sensors in the sensor description. the unused components need to
  % be removed in a second step
  sens = ft_apply_montage(data.(sensfield), montage, 'keepunused', 'yes', 'balancename', 'invcomp', 'feedback', cfg.feedback);
  
  % there could have been sequential subspace projections, so the
  % invcomp-field may have been renamed into invcompX. If this it the case,
  % take the one with the highest suffix
  invcompfield = 'invcomp';
  if ~isfield(sens.balance, 'invcomp')
    for k = 10:-1:1
      if isfield(sens.balance, ['invcomp',num2str(k)])
        invcompfield = [invcompfield,num2str(k)];
        break;
      end
    end
  end
  
  % remove the unused channels from the grad/elec
  [junk, remove]    = match_str(comp.label, sens.label);
  sens.tra(remove,:) = [];
  sens.label(remove) = [];
  sens.chanpos(remove,:) = [];
  if isfield(sens, 'chanori')
    sens.chanori(remove,:) = [];
  end
  
  % remove the unused components from the balancing and from the tra
  [junk, remove]    = match_str(comp.label, sens.balance.(invcompfield).labelnew);
  sens.balance.(invcompfield).tra(remove, :)   = [];
  sens.balance.(invcompfield).labelnew(remove) = [];
  data.(sensfield)  = sens;
end

if istlck
  % convert the raw structure back into a timelock structure
  data = ft_checkdata(data, 'datatype', 'timelock');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   comp data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
