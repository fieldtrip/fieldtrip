% FT_PREAMBLE_LOADVAR is a helper script that optionally loads one or
% multiple fieldtrip data structures from mat files on disk, as an
% alternative to the user specifying the data structures as input variables
% to the calling function. This makes use of the cfg.inputfile variable.
%
% Furthermore, this function computes the MD5 hash of the input data structures for
% provenance.
%
% Use as
%   ft_preamble loadvar data
%   ft_preamble loadvar source mri
%   ft_preamble loadvar varargin
%
% See also FT_POSTAMBLE_SAVEVAR, FT_PREAMBLE_PROVENANCE

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
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

% the name of the variables are passed in the preamble field
global ft_default

% use an anonymous function
assign = @(var, val) assignin('caller', var, val);

if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)
  
  % the input data should be read from file
  if (nargin>1)
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
      mutexlock(cfg.inputlock);
    end
    
    if isequal(ft_default.preamble, {'varargin'}) && ~iscell(cfg.inputfile)
      % this should be a cell-array, oterwise it cannot be assigned to varargin
      cfg.inputfile = {cfg.inputfile};
    end
    
    if iscell(cfg.inputfile)
      if isequal(ft_default.preamble, {'varargin'})
        % read multiple inputs and copy them into varargin
        tmp = {};
        for i=1:length(cfg.inputfile)
          tmp{i} = loadvar(cfg.inputfile{i}, 'data');
        end % for
        assign('varargin', tmp);
        clear i tmp
      else
        % ft_default.preamble is a cell-array containing the variable names
        for i=1:length(cfg.inputfile)
          assign(ft_default.preamble{i}, loadvar(cfg.inputfile{i}, ft_default.preamble{i}));
        end % for
        clear i
      end
    else
      % ft_default.preamble{1} contains the variable name
      assign(ft_default.preamble{1}, loadvar(cfg.inputfile, ft_default.preamble{1}));
    end
    
    if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
      mutexunlock(cfg.inputlock);
    end
  end
  
end % if inputfile

% the following section deals with tracking the information about the input data structures
% the corresponding section for the output data structures is in ft_postamble_history

if isfield(cfg, 'trackdatainfo') && istrue(cfg.trackdatainfo)
  if isequal(ft_default.preamble, {'varargin'})
    for i=1:length(varargin)
      % store the hash for each input argument
      cfg.datainfo.input{i} = hashvar(varargin{i});
    end
  else
    % ft_default.preamble is a cell-array containing the variable names
    for i=1:length(ft_default.preamble)
      if exist(ft_default.preamble{i}, 'var')
        cfg.datainfo.input{i} = eval(sprintf('hashvar(%s)', ft_default.preamble{i}));
      else
        % the function can have optional inputs that are unspecified, e.g. data for ft_preprocessing
        cfg.datainfo.input{i} = [];
      end
    end
  end
end
