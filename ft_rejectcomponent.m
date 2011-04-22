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

ft_defaults

% set defaults
if ~isfield(cfg, 'component'), cfg.component = [];      end
if ~isfield(cfg, 'inputfile'),    cfg.inputfile = [];           end
if ~isfield(cfg, 'outputfile'),   cfg.outputfile = [];          end

if nargin==3 
  %ntrials = length(data.trial);
  data    = ft_checkdata(data, 'datatype', 'raw');
  label   = data.label;
  hasdata = 1;
elseif nargin==2 
  %ntrials = length(comp.trial);
  label   = comp.topolabel;
  hasdata = 0;
elseif nargin<2 % only cfg is given; inputfile is expected
  comp = loadvar(cfg.inputfile, 'comp');
  %ntrials = length(comp.trial);
  label   = comp.topolabel;
  hasdata = 0;
end

comp    = ft_checkdata(comp, 'datatype', 'comp');
ncomps  = length(comp.label);

if min(cfg.component)<1
  error('you cannot remove components that are not present in the data');
end

if max(cfg.component)>ncomps
  error('you cannot remove components that are not present in the data');
end

% set the rejected component amplitudes to zero 
fprintf('removing %d components\n', length(cfg.component)); 
fprintf('keeping %d components\n',  ncomps-length(cfg.component));

%create a projection matrix by subtracting the subspace spanned by the 
%topographies of the to-be-removed components from identity
[seldat, selcomp] = match_str(label, comp.topolabel);

if length(seldat)~=length(label) && nargin==3,
  warning('the subspace projection is not guaranteed to be correct for non-orthogonal components');
end

if hasdata,
  topo     = comp.topo(selcomp,:);
  invtopo  = pinv(topo);
  tra      = eye(length(selcomp)) - topo(:, cfg.component)*invtopo(cfg.component, :);
  %I am not sure about this, but it gives comparable results to the ~hasdata case
  %when comp contains non-orthogonal (=ica) topographies, and contains a complete decomposition
  
  %we are going from data to components, and back again
  labelorg = comp.topolabel(selcomp);
  labelnew = comp.topolabel(selcomp);

  keepunused = 'yes'; %keep the original data which are not present in the mixing provided
  
else
  topo = comp.topo(selcomp, :);
  topo(:, cfg.component) = 0;
  tra      = topo;
  
  %we are going from components to data
  labelorg = comp.label;
  labelnew = comp.topolabel(selcomp);
  
  %create data structure
  if hasdata && isfield(data, 'trialinfo'),  trialinfo  = data.trialinfo;  end
  if hasdata && isfield(data, 'sampleinfo'), sampleinfo = data.sampleinfo; end 
  data         = [];
  data.trial   = comp.trial;
  data.time    = comp.time;
  data.label   = comp.label;
  data.fsample = comp.fsample;
  if isfield(comp, 'grad'), data.grad       = comp.grad;  end
  if isfield(comp, 'elec'), data.elec       = comp.elec;  end
  if exist('trialinfo',  'var'),   data.trialinfo  = trialinfo;  end
  if exist('sampleinfo', 'var'),   data.sampleinfo = sampleinfo; end
  
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
data             = ft_apply_montage(data, montage, 'keepunused', keepunused);
if isfield(data, 'grad') || (isfield(data, 'elec') && isfield(data.elec, 'tra')),
  if isfield(data, 'grad')
    sensfield = 'grad';
  else
    sensfield = 'elec';
  end
  data.(sensfield) = ft_apply_montage(data.(sensfield), montage, 'keepunused', keepunused, 'balancename', 'invcomp');
else
  %warning('the gradiometer description does not match the data anymore');
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add the version details of this function call to the configuration 
cfg.version.name = mfilename('fullpath'); 
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.version.matlab = version();

if ~hasdata 
  % remember the configuration details of the input data
  if isfield(comp, 'cfg'), cfg.previous = comp.cfg; end
  % copy the sampleinfo into the output
  if isfield(comp, 'sampleinfo')
    data.sampleinfo = comp.sampleinfo;
  end
  % copy the trialinfo into the output
  if isfield(comp, 'trialinfo')
    data.trialinfo = comp.trialinfo;
  end
elseif hasdata
  if isfield(comp, 'cfg'), cfg.previous{1} = comp.cfg; end
  if isfield(comp, 'cfg'), cfg.previous{2} = data.cfg; end
end

% keep the configuration in the output
data.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', data); % use the variable name "data" in the output file
end

