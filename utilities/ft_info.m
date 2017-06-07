function ft_info (feedback, varargin)
% FT_INFO will print an info messages to the command window unsing fprintf, 
% if feedback is enabled. Otherwise the message will be discarded.
%
% Use as
%   FT_INFO(feedback, FORMAT, A, ...) applies the FORMAT to all elements of 
%   array A and any additional array arguments in column order, and
%   displays the result on the screen, if feedback is set to 'yes'.
%
% The input data can contain the following values
%   feedback      = 'yes' or 'no'
%
% SEE also FPRINTF

% Copyright (C) 2017, Robert Oostenveld, Daniel Matthes
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


if strcmp(feedback, 'yes')
  fprintf('info: ');
  fprintf(varargin{:});
else
  % do nothing
end