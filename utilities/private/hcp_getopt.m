function [opts, args] = hcp_getopt(varargin)

% HCP_GETOPT parses the command line to the megconnectome executable
% application, separating the options that start with -- from the file
% names of the scripts to be executed.
%
% Use as
%   megconnectome.exe --option1 arg1 --option2 arg2 scriptA.m scriptB.m
% splits the command line arguments into a cell array with key-value pairs
% and a cell array with the filenames.
%
% In this example the hcp_getopt function returns
%   opts = {'option1', arg1, 'option2', arg2};
%   args = {'scriptA.m', 'scriptB.m'}
%
% See also FT_GETOPT

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
%
% This file is part of megconnectome.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.

opts = {};
args = {};

i = 1;
while (i<=nargin)
  if ischar(varargin{i}) && strncmp('--', varargin{i}, 2)
    % assume that it is a --option
    opts{end+1} = varargin{i}(3:end); % remove the -- from the option
    try
      opts{end+1} = varargin{i+1};      % the value
    catch
      opts{end+1} = [];                 % there are not enough values
    end
    i = i + 2;
  elseif ischar(varargin{i}) && strncmp('-', varargin{i}, 1)
    % assume that it is a -option
    opts{end+1} = varargin{i}(2:end); % remove the - from the option
    try
      opts{end+1} = varargin{i+1};      % the value
    catch
      opts{end+1} = [];                 % there are not enough values
    end
    i = i + 2;
  else
    % assume that it is a file name
    args{end+1} = varargin{i};
    i = i + 1;
  end
end

