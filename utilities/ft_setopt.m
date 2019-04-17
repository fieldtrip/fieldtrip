function opt = ft_setopt(opt, key, val)

% FT_SETOPT assigns a value to an configuration structure or to a cell-array
% with key-value pairs. It will overwrite the option if already present, or
% append the option if not present.
%
% Use as
%   s = ft_setopt(s, key, val)
% where s is a structure or a cell-array.
%
% See also FT_GETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
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

if isa(opt, 'struct') || isa(opt, 'config')

  % just replace or add the option
  opt.(key) = val;

elseif isa(opt, 'cell')
  % determine whether the key already exists
  fn = opt(1:2:end);
  sel = find(strcmp(key, fn));
  if isempty(sel)
    % append it
    opt{end+1} = key;
    opt{end+1} = val;
  elseif length(sel)==1
    % replace the current value
    keyindex = 2*sel-1;
    valindex = keyindex+1;
    opt{valindex} = val;
  elseif length(sel)>1
    % first remove all occurences
    keyindex = 2*sel-1;
    valindex = keyindex+1;
    opt([keyindex valindex]) = [];
    % then append it
    opt{end+1} = key;
    opt{end+1} = val;
  end
end % isstruct or iscell
