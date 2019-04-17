function [varargout] = webwrite(varargin)

% WEBWRITE is a drop-in replacement for the function with the same
% name that was introduced in MATLAB 2014b. This function is only
% partially compatible with the original.
%
% Usage:
%   status = webwrite(url, [data], [options])
%
% This function posts data (which can be a struct or a string) to url.
% The options argument is currently ignored.
%
% Unlike builtin MATLAB webwrite, this function
% - sets headers always to 'Content-Type: application/json; charset=UTF-8'
%   irrespective of the contents of the options argument.
% - the status returned is zero if succesful, and nonzero otherwise
%
% This requires that curl is available on the command-line.
% If curl is not available, or an other error occurs, then this function
% raises an error.

% Copyright (C) 2017-2019, Robert Oostenveld, Nikolaas N. Oosterhof
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see https://github.com/fieldtrip/fieldtrip/issues/899

if exist(mfilename, 'builtin') || any(strncmp(which(mfilename, '-all'), matlabroot, length(matlabroot)) & cellfun(@isempty, regexp(which(mfilename, '-all'), fullfile('private', mfilename))))
  % remove this directory from the path
  p = fileparts(mfilename('fullpath'));
  warning('removing %s from your path, see http://bit.ly/2SPPjUS', p);
  rmpath(p);
  % call the original MATLAB function
  if exist(mfilename, 'builtin')
    [varargout{1:nargout}] = builtin(mfilename, varargin{:});
  else
    [varargout{1:nargout}] = feval(mfilename, varargin{:});
  end
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is where the actual replacement code starts
% function status = webwrite(url, data, options)

% deal with the input arguments
if nargin==1
  [url               ] = deal(varargin{1:1});
elseif nargin==2
  [url, data         ] = deal(varargin{1:2});
elseif nargin==3
  [url, data, options] = deal(varargin{1:3});
else
  error('incorrect number of input arguments')
end

if nargin<2
  data = '';
end

if isstruct(data)
  % convert it to JSON string
  data = struct2json(data);
  assert(ischar(data));

  to_remove=sprintf('[\n ]'); % newlines and spaces
  data = regexprep(data, to_remove, '');
end

if ~ischar(data)
  error('unsupported data: input must be string or struct');
end

cmd = sprintf(['curl '...
                '-s '...                                % quiet (no output)
                '-H "Content-Type: application/json;'...  % JSON data
                    'charset=UTF-8" '...                % charset
                '-X POST '...                           % method
                '-d ''%s'' %s'], ...                    % data and URL
                data, url);

[status, str] = system(cmd);

if status
  error('An error occured when running this command: ''%s'':\n%s',cmd,str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deal with the output arguments

varargout = {status};
