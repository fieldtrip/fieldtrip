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
ft_preamble loadvar comp data
ft_preamble provenance comp data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.component       = ft_getopt(cfg, 'component',  []);
cfg.demean          = ft_getopt(cfg, 'demean',     'yes');
cfg.feedback        = ft_getopt(cfg, 'feedback',   'text');
cfg.updatesens      = ft_getopt(cfg, 'updatesens', 'yes');

% the data can be passed as input arguments or can be read from disk
hascomp = exist('comp', 'var');
hasdata = exist('data', 'var');

if hascomp && hasdata
  % check if the input data is valid for this function
  istlck  = ft_datatype(data, 'timelock');  % this will be temporary converted into raw
  data    = ft_checkdata(data, 'datatype', 'raw');
  comp    = ft_checkdata(comp, 'datatype', 'comp');
  label   = data.label;
  nchans  = length(data.label);
  ncomps  = length(comp.label);
elseif hascomp
  % check if the input data is valid for this function
  istlck  = ft_datatype(comp, 'timelock');  % this will be temporary converted into raw
  comp    = ft_checkdata(comp, 'datatype', 'raw+comp');
  label   = comp.topolabel;
  ncomps  = length(comp.label);
else
  ft_error('incorrect number of input arguments');
end

% cfg.component can be indicated by number or by label
cfg.component = ft_channelselection(cfg.component, comp.label);
reject = match_str(comp.label, cfg.component);

if isempty(reject)
  ft_warning('no components were selected for rejection');
end

if min(reject)<1
  ft_error('you cannot remove components that are not present in the data');
end

if max(reject)>ncomps
  ft_error('you cannot remove components that are not present in the data');
end

if hasdata && strcmp(cfg.demean, 'yes')
  % optionally perform baseline correction on each trial
  ft_info('baseline correcting data \n');
  for trial=1:numel(data.trial)
    data.trial{trial} = ft_preproc_baselinecorrect(data.trial{trial});
  end
end

% set the rejected component amplitudes to zero
ft_info('removing %d components\n', length(reject));
ft_info('keeping %d components\n',  ncomps-length(reject));

% create a projection matrix by subtracting the subspace spanned by the
% topographies of the to-be-removed components from identity
[seldat, selcomp] = match_str(label, comp.topolabel);

if hasdata && length(seldat)~=length(label)
  ft_warning('the subspace projection is not guaranteed to be correct for non-orthogonal components');
end

if hasdata
  mixing   = comp.topo(selcomp,:);
  unmixing = comp.unmixing(:,selcomp);
  
  % I am not sure about this, but it gives comparable results to the ~hasdata case
  % when comp contains non-orthogonal (=ica) topographies, and contains a complete decomposition
  
  montage     = [];
  montage.tra = eye(length(selcomp)) - mixing(:, reject)*unmixing(reject, :);
  % we are going from data to components, and back again
  montage.labelold = comp.topolabel(selcomp);
  montage.labelnew = comp.topolabel(selcomp);
  
  keepunused = 'yes'; % keep the original data which are not present in the mixing provided
  bname = 'reject';
  
else
  mixing = comp.topo(selcomp, :);
  mixing(:, reject) = 0;
  
  montage     = [];
  montage.tra = mixing;
  % we are going from components to data
  montage.labelold = comp.label;
  montage.labelnew = comp.topolabel(selcomp);
  
  keepunused = 'no'; % don't need to keep the original rejected components
  bname = 'invcomp';
  
  % create the initial data structure, remove all component details
  data = keepfields(comp, {'trial', 'time', 'label', 'fsample', 'grad', 'elec', 'opto', 'trialinfo', 'sampleinfo'});
end % if hasdata

% apply the linear projection to the data
data = ft_apply_montage(data, montage, 'keepunused', keepunused, 'feedback', cfg.feedback, 'showcallinfo', cfg.showcallinfo);

sensfield = cell(0,1);
if isfield(data, 'grad')
  sensfield{end+1} = 'grad';
end
if isfield(data, 'elec')
  sensfield{end+1} = 'elec';
end
if isfield(data, 'opto')
  sensfield{end+1} = 'opto';
end

% apply the linear projection also to the sensor description
if ~isempty(sensfield)
  if  strcmp(cfg.updatesens, 'yes')

    for m = 1:numel(sensfield)
      ft_info('also applying the backprojection matrix to the %s structure\n', sensfield{m});
      
      % the balance field is needed to keep the sequence of linear projections
      if ~isfield(data.(sensfield{m}), 'balance')
        data.(sensfield{m}).balance.current = 'none';
      end
      
      % keepunused = 'yes' is required to get back e.g. reference or otherwise
      % unused sensors in the sensor description. the unused components need to
      % be removed in a second step
      sens = ft_apply_montage(data.(sensfield{m}), montage, 'keepunused', 'yes', 'balancename', bname, 'feedback', cfg.feedback);
      
      % remove the unused channels from the grad/elec/opto
      [junk, remove]    = match_str(comp.label, sens.label);
      sens.tra(remove,:) = [];
      sens.label(remove) = [];
      sens.chanpos(remove,:) = [];
      if isfield(sens, 'chanori')
        sens.chanori(remove,:) = [];
      end
      
      % there could have been sequential subspace projections, so the
      % invcomp-field may have been renamed into invcompX. If this it the case,
      % take the one with the highest suffix
      invcompfield = bname;
      if  ~isfield(sens.balance, invcompfield)
        for k = 10:-1:1
          if isfield(sens.balance, [bname num2str(k)])
            invcompfield = [invcompfield num2str(k)];
            break;
          end

        end
      end
      
      % remove the unused components from the balancing
      [junk, remove]    = match_str(comp.label, sens.balance.(invcompfield).labelnew);
      sens.balance.(invcompfield).tra(remove, :)   = [];
      sens.balance.(invcompfield).labelnew(remove) = [];
      data.(sensfield{m})  = sens;
    end
    
  else
    for m = 1:numel(sensfield)
      ft_info('not applying the backprojection matrix to the %s structure\n', sensfield{m});
      % simply copy it over
      comp.(sensfield{m}) = data.(sensfield{m});
    end
  end
end % if sensfield

if istlck
  % convert the raw structure back into a timelock structure
  data = ft_checkdata(data, 'datatype', 'timelock');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   comp data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data
