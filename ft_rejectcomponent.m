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
% An optional input argument data can be provided. In that case
% componentanalysis will do a subspace projection of the input data
% onto the space which is spanned by the topographies in the unmixing
% matrix in comp, after removal of the artifact components.
%
% The configuration should contain
%   cfg.component = list of components to remove, e.g. [1 4 7]
%   cfg.demean    = 'no' or 'yes', whether to demean the input data (default = 'yes')
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

% Copyright (C) 2005-2009, Robert Oostenveld
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
ft_preamble provenance
ft_preamble trackconfig
ft_preamble debug
ft_preamble loadvar comp data

% set defaults
cfg.component  = ft_getopt(cfg, 'component',  []);
cfg.demean     = ft_getopt(cfg, 'demean',    'yes');
cfg.feedback   = ft_getopt(cfg, 'feedback',  'text');

% the data can be passed as input arguments or can be read from disk
nargin = 1;
nargin = nargin + exist('comp', 'var');
nargin = nargin + exist('data', 'var');

if nargin==3
  data    = ft_checkdata(data, 'datatype', 'raw');
  label   = data.label;
  hasdata = 1;
  nchans  = length(data.label);
elseif nargin==2
  label   = comp.topolabel;
  hasdata = 0;
else
  error('incorrect number of input arguments');
end

comp    = ft_checkdata(comp, 'datatype', 'comp');
ncomps  = length(comp.label);


if min(cfg.component)<1
  error('you cannot remove components that are not present in the data');
end

if max(cfg.component)>ncomps
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
fprintf('removing %d components\n', length(cfg.component));
if ~hasdata,
  fprintf('keeping %d components\n',  ncomps-length(cfg.component));
else
  fprintf('keeping %d components\n',  nchans-length(cfg.component));
end

%create a projection matrix by subtracting the subspace spanned by the
%topographies of the to-be-removed components from identity
[seldat, selcomp] = match_str(label, comp.topolabel);

if length(seldat)~=length(label) && nargin==3,
  warning('the subspace projection is not guaranteed to be correct for non-orthogonal components');
end

if hasdata
  mixing = comp.topo(selcomp,:);
  unmixing = comp.unmixing(:,selcomp);
  tra = eye(length(selcomp)) - mixing(:, cfg.component)*unmixing(cfg.component, :);
  %I am not sure about this, but it gives comparable results to the ~hasdata case
  %when comp contains non-orthogonal (=ica) topographies, and contains a complete decomposition
  
  %we are going from data to components, and back again
  labelorg = comp.topolabel(selcomp);
  labelnew = comp.topolabel(selcomp);
  
  keepunused = 'yes'; %keep the original data which are not present in the mixing provided
  
else
  mixing = comp.topo(selcomp, :);
  mixing(:, cfg.component) = 0;
  tra = mixing;
  
  %we are going from components to data
  labelorg = comp.label;
  labelnew = comp.topolabel(selcomp);
  
  %create data structure
  data         = [];
  data.trial   = comp.trial;
  data.time    = comp.time;
  data.label   = comp.label;
  data.fsample = comp.fsample;
  if isfield(comp, 'grad'), data.grad       = comp.grad;  end
  if isfield(comp, 'elec'), data.elec       = comp.elec;  end
  if isfield(comp, 'trialinfo'),   data.trialinfo  = comp.trialinfo;  end
  if isfield(comp, 'sampleinfo'),   data.sampleinfo = comp.sampleinfo; end
  
  keepunused = 'no'; %don't need to keep the original rejected components
end

%OLD CODE
% recontruct the trials
%for i=1:ntrials
%  data.trial{i} = projector * data.trial{i}(seldat,:);
%end
%data.label = data.label(seldat);

%create montage and apply this to data and grad
montage          = [];
montage.tra      = tra;
montage.labelorg = labelorg;
montage.labelnew = labelnew;
data             = ft_apply_montage(data, montage, 'keepunused', keepunused, 'feedback', cfg.feedback);

if isfield(data, 'grad') || (isfield(data, 'elec') && isfield(data.elec, 'tra')),
  if isfield(data, 'grad')
    sensfield = 'grad';
  else
    sensfield = 'elec';
  end
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
  %data.(sensfield)  = ft_apply_montage(data.(sensfield), montage, 'keepunused', 'no', 'balancename', 'invcomp');
else
  %warning('the gradiometer description does not match the data anymore');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance
if nargin==2
  ft_postamble previous comp
elseif nargin==3
  ft_postamble previous comp data
end
ft_postamble history data
ft_postamble savevar data
