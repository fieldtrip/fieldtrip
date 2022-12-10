function input = handle_edit_input(text)

% HANDLE_EDIT_INPUT deals with user-entered input in the GUI. This is used to select
% channels and/or trials in FT_REJECTVISUAL and to select channels in FT_DATABROWSER
%
% The input text can consist of a string such as
%   1 2 3 4
%   1:4
%   [1 2 3 4]
%   [1:4]
% This is converted in a list of numbers.
%
% The input text can also consist of a single non-numeric string or a string that
% represents a cell-array of strings such as
%   all
%   {'MEG', '-MR*'}

% Copyright (C) 2022, Robert Oostenveld
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

if isempty(text)
  input = [];
elseif all(ismember(text, '0123456789 \t'))
  input = regexp(text, '\s+', 'split');
  input = cellfun(@str2num, input);
elseif ismember(text(1), '0123456789')
  input = eval(text); % it could be a single number or a string like 1:10
elseif startsWith(text, '[') && endsWith(text, ']')
  input = eval(text);
elseif startsWith(text, '{') && endsWith(text, '}')
  input = eval(text); % this allows for {'MEG', '-MR*'}
else
  input = text;
end
