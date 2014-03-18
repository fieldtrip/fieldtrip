function data = ft_deidentifydata(cfg, data)

% FT_DEIDENTIFYDATA clears the value of potentially identifying fields in
% the provenance information, i.e., it updates the configuration structure
% and history that is maintained by FieldTrip in the cfg field.

% Use as
%  outdata = ft_deidentifydata(cfg, indata)
% where indata is any FieldTrip data structure and cfg is a configuration
% structure that should contain
%
%  cfg.keep      = cell-array with strings, fields to keep (default = {})
%  cfg.keep      = cell-array with strings, fields to remove (default = {})
%  cfg.uncertain = string, what to do with unknown fields, can be 'ask',
%                  'keep', 'remove' (default = 'ask')
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_ANALYSISPIPELINE

% Copyright (C) 2014, Robert Oostenveld, DCCN
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
ft_defaults                 % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init            % this will reset warning_once and show the function help if nargin==0 and return an error
ft_preamble provenance      % this records the time and memory usage at teh beginning of the function
ft_preamble trackconfig     % this converts the cfg structure in a config object, which tracks the cfg options that are being used
ft_preamble debug           % this allows for displaying or saving the function name and input arguments upon an error
ft_preamble loadvar data    % this reads the input data in case the user specified the cfg.inputfile option

% get the options
cfg.keep      = ft_getopt(cfg, 'keep',      {});
cfg.remove    = ft_getopt(cfg, 'remove',    {});
cfg.uncertain = ft_getopt(cfg, 'uncertain', 'ask');


% process the data using a recursive helper function
if ~isfield(data, 'cfg')
  warning('the input data has no provenance information, nothing to do');
else
  data.cfg = struct(data.cfg); % ensure that it is a structure
  
  % these should be column arrays
  originalkeep = cfg.keep(:);
  originalremove = cfg.remove(:);
  
  % these two variables are recursively updated by the subfunction
  actualkeep   = {};
  actualremove = {};
  
  [data.cfg, actualkeep, actualremove] = parsestructure(data.cfg, cfg.keep, cfg.remove, actualkeep, actualremove, cfg.uncertain);
  
  % update the output configuration with the actual fields
  cfg.keep   = actualkeep;
  cfg.remove = actualremove;
  
end

% deal with the output
ft_postamble debug            % this clears the onCleanup function used for debugging in case of an error
ft_postamble trackconfig      % this converts the config object back into a struct and can report on the unused fields
ft_postamble provenance       % this records the time and memory at the end of the function, prints them on screen and adds this information together with the function name and matlab version etc. to the output cfg
ft_postamble previous data    % this copies the datain.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history data     % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg
ft_postamble savevar data     % this saves the output data structure to disk in case the user specified the cfg.outputfile option

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s, actualkeep, actualremove, uncertain] = parsestructure(s, originalkeep, originalremove, actualkeep, actualremove, uncertain)
fn = fieldnames(s);

% concatenate the known fields to keep and remove
originalkeep   = unique(cat(1, originalkeep, actualkeep));
originalremove = unique(cat(1, originalremove, actualremove));

for i=1:length(fn)
  
  if isstruct(s.(fn{i}))
    % recurse into all sub-structures
    [s.(fn{i}), actualkeep, actualremove, uncertain] = parsestructure(s.(fn{i}), originalkeep, originalremove, actualkeep, actualremove, uncertain);
    
  elseif ~ischar(s.(fn{i}))
    % fprintf('keeping field %s with non-string value\n', fn{i});
    
  elseif any(strcmp(fn{i}, originalkeep)) || strcmp(uncertain, 'keep')
    fprintf('keeping %s = ''%s''\n', fn{i},  s.(fn{i}));
    actualkeep = cat(1, actualkeep, fn{i});
    
  elseif any(strcmp(fn{i}, originalremove)) || strcmp(uncertain, 'remove')
    fprintf('removing %s = ''%s''\n', fn{i},  s.(fn{i}));
    s.(fn{i}) = 'removed by ft_deidentifydata';
    actualremove = cat(1, actualremove, fn{i});
    
  elseif strcmp(uncertain, 'ask')
    fprintf('\nUncertain about %s = ''%s''\n', fn{i},  s.(fn{i}));
    question = 'Should this field be (k)ept, (r)emoved, [K]eep all, [R]emove all? ';
    
    answer = 'x';
    while ~ismember(answer, {'k', 'r', 'K', 'R'})
      answer = smartinput(question, answer);
    end
    
    switch answer
      case 'k'
        fprintf('keeping %s = ''%s''\n', fn{i},  s.(fn{i}));
        actualkeep = cat(1, actualkeep, fn{i});
      case 'r'
        fprintf('removing %s = ''%s''\n', fn{i},  s.(fn{i}));
        actualremove = cat(1, actualremove, fn{i});
        s.(fn{i}) = 'removed by ft_deidentifydata';
      case 'K'
        fprintf('keeping %s = ''%s''\n', fn{i},  s.(fn{i}));
        actualkeep = cat(1, actualkeep, fn{i});
        uncertain = 'keep'; % for the next time
      case 'R'
        fprintf('removing %s = ''%s''\n', fn{i},  s.(fn{i}));
        actualremove = cat(1, actualremove, fn{i});
        s.(fn{i}) = 'deidentified';
        uncertain = 'remove'; % for the next time
        
    end % switch
    
  else
    error('Cannot decide what to do with the field %s', fn{i});
  end
end % for

% ensure that each field is only present once
actualkeep   = unique(actualkeep);
actualremove = unique(actualremove);

