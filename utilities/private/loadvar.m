function value = loadvar(filename, varname)

% LOADVAR is a helper function for cfg.inputfile
%
% See also SAVEVAR

% Copyright (C) 2010, Robert Oostenveld
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

assert(ischar(filename), 'file name should be a string');

if nargin<2
  ft_info('reading variable from file ''%s''\n', filename);
else
  assert(ischar(varname), 'variable name should be a string');
  ft_info('reading ''%s'' from file ''%s''\n', varname, filename);
end

% note that this sometimes fails, returning an empty var
% this is probably due to MATLAB filename and MATLAB version issues
var = whos('-file', filename);

if isempty(var) && nargin==1
  filecontent = load(filename); % read everything from the file, regardless of how the variables are called
  varname = fieldnames(filecontent);
  if length(varname)==1
    % the one variable in the file will be returned
    value = filecontent.(varname{1});
    clear filecontent
  else
    ft_error('cannot read an unspecified variable in case of a file containing multiple variables');
  end

elseif length(var)==1
  filecontent = load(filename); % read the one variable in the file, regardless of how it is called
  value       = filecontent.(var.name);
  clear filecontent

else
  filecontent = load(filename, varname);
  value       = filecontent.(varname);  % read the variable named according to the input specification
  clear filecontent
end

