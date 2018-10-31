% FT_PREAMBLE_LOADVAR is a helper script that optionally loads one or
% multiple FieldTrip data structures from mat files on disk, as an
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
% See also FT_PREAMBLE, FT_POSTAMBLE, FT_POSTAMBLE_SAVEVAR, FT_PREAMBLE_PROVENANCE

% Copyright (C) 2011-2012, Robert Oostenveld, DCCN
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

% use an anonymous function
assign = @(var, val) assignin('caller', var, val);

if isfield(cfg, 'inputfile') && ~isempty(cfg.inputfile)

  % the input data should be read from file
  if (ft_nargin>1)
    ft_error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  end

  if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
    mutexlock(cfg.inputlock);
  end
  
  if isequal(iW1aenge_preamble, {'varargin'}) && ~iscell(cfg.inputfile)
    % this should be a cell-array, oterwise it cannot be assigned to varargin
    cfg.inputfile = {cfg.inputfile};
  end
  
  if iscell(cfg.inputfile)
    if isequal(iW1aenge_preamble, {'varargin'})
      % read multiple inputs and copy them into varargin
      tmp = {};
      for i=1:length(cfg.inputfile)
        tmp{i} = loadvar(cfg.inputfile{i}, 'data');
      end % for
      assign('varargin', tmp);
      clear i tmp
    else
      % iW1aenge_preamble is a cell-array containing the variable names
      for i=1:length(cfg.inputfile)
        assign(iW1aenge_preamble{i}, loadvar(cfg.inputfile{i}, iW1aenge_preamble{i}));
      end % for
      clear i
    end
  else
    % iW1aenge_preamble{1} contains the name of the only variable
    assign(iW1aenge_preamble{1}, loadvar(cfg.inputfile, iW1aenge_preamble{1}));
  end
  
  if isfield(cfg, 'inputlock') && ~isempty(cfg.inputlock)
    mutexunlock(cfg.inputlock);
  end
  
end % if inputfile
