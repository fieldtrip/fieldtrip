function dimsiz = getdimsiz(data, field, numdim)

% GETDIMSIZ
%
% Use as
%   dimsiz = getdimsiz(data, field)
% or
%   dimsiz = getdimsiz(data, field, numdim)
%
% MATLAB will not return the size of a  field in the data structure that has trailing
% singleton dimensions, since those are automatically squeezed out. With the optional
% numdim parameter you can specify how many dimensions the data element has. This
% will result in the trailing singleton dimensions being added to the output vector.
%
% Example use
%   dimord = getdimord(datastructure, fieldname);
%   dimtok = tokenize(dimord, '_');
%   dimsiz = getdimsiz(datastructure, fieldname, numel(dimtok));
%
% See also GETDIMORD, GETDATFIELD

% Copyright (C) 2014-2019, Robert Oostenveld
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

% deal with trailing singleton dimensions
if nargin<3
  numdim = 0;
end

% deal with data located in a nested sub-structure
if ~isfield(data, field) && isfield(data, 'avg') && isfield(data.avg, field)
  field = ['avg.' field];
elseif ~isfield(data, field) && isfield(data, 'trial') && isfield(data.trial, field)
  field = ['trial.' field];
elseif ~isfield(data, field)
  ft_error('field "%s" not present in data', field);
end

if strncmp(field, 'avg.', 4)
  prefix = [];
  field = field(5:end); % strip the avg
  data.(field) = data.avg.(field); % move the avg into the main structure
  data = rmfield(data, 'avg');
elseif strncmp(field, 'trial.', 6)
  prefix = numel(data.trial);
  field = field(7:end); % strip the trial
  data.(field) = data.trial(1).(field); % move the first trial into the main structure
  data = rmfield(data, 'trial');
else
  prefix = [];
end

dimsiz = cellmatsize(data.(field));

% add nrpt in case of source.trial
dimsiz = [prefix dimsiz];

% deal with trailing singleton dimensions
dimsiz(end+1:numdim) = 1;

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to determine the size of data representations like {pos}_ori_time
% FIXME this will fail for {xxx_yyy}_zzz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function siz = cellmatsize(x)
if iscell(x)
  if isempty(x)
    siz = 0;
    return % nothing else to do
  elseif isvector(x)
    cellsize = numel(x);          % the number of elements in the cell-array
  else
    cellsize = size(x);
    x = x(:); % convert to vector for further size detection
  end
  [dum, indx] = max(cellfun(@numel, x));
  matsize = size(x{indx});      % the size of the content of the cell-array
  siz = [cellsize matsize];     % concatenate the two
else
  siz = size(x);
end
end % function cellmatsize
