function [data] = swapmemfile(data)

% SWAPMEMFILE swaps a variable from file into memory and clears it 
% again from the memory on the subsequent call
%
% Use with extreme caution!

% This function depends on PREPARE_TIMEFREQ_DATA which has the following options:
% cfg.avgoverchan
% cfg.avgoverfreq
% cfg.avgovertime
% cfg.channel            
% cfg.channelcmb         
% cfg.datarepresentation
% cfg.frequency
% cfg.latency            
% cfg.precision
% cfg.previous
% cfg.version

% Copyright (C) 2004, Robert Oostenveld
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

% this variable will be empty at the first call
persistent file

if nargin<1 || isempty(data)
  % (re)initialize by clearing the persistent variable
  file = [];
  data = [];
elseif ~isempty(file)
  % data is present in memory and originates from a file
  % remove data from memory and replace it by the filename
  fprintf('swapping %s out of memory\n', file);
  data = file;
  file = [];
elseif isempty(file) && (~ischar(data) || ~filetype(data, 'matlab'))
  % data is already present in memory, and does not originate from a file
  % do nothing
elseif isempty(file) && ischar(data) && ft_filetype(data, 'matlab')
  % data is not present in memory but should be read from a MATLAB file
  % remember the filename from which the data originates
  file = data;
  fprintf('swapping %s into memory\n', file);
  data = load(file);
  % select the variable of interest, discard the rest
  str = fieldnames(data);
  if length(str)==1
    % select the single variable that is inside the MATLAB file
    data = getfield(data, str{1});
  else
    error('Matlab file should contain only one variable');
  end
else
  error('Unknown combination of input variable and file history');
end
