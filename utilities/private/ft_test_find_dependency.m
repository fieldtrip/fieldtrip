function [outlist, depmat] = find_dependency(inlist, varargin)

% FIND_DEPENDENCY creates a dependency matrix for direct dependencies of test scripts
% 
% Use as
%    [outlist, depmat] = find_dependency(inlist)
% where 
%   inlist   = Nx1 cell-array, descibes the rows and lists the test scripts
%   outlist  = 1xM cell-array, describes the columns and lists the
%   dependencies
%   depmat   = NxM matrix, see below
% 
% The depmat contains the following values:
%  - 0 if there is no dependency
%  - 2 for a direct dependency
% 
% See also UPDATE_DEPENDENCY, UNTESTED_FUNCTIONS

% Copyright (C) 2023, Konstantinos Tsilimparis
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

for i=1:length(inlist)  
  [p, f, x] = fileparts(inlist{i});
  inlist{i} = f; % Remove path and file extension

  try
    fprintf('processing "%s"\n', inlist{i});
    dep{i} = matlab.codetools.requiredFilesAndProducts(inlist{i}, 'toponly'); % Determine direct dependencies

    indices           = ~contains(dep{i}, inlist{i});    % Remove self-dependencies
    indices_dccnpath  = ~contains(dep{i}, "dccnpath");   % Remove the dependency on dccnpath, since dccnpath is widely used
    indices_ft_getopt = ~contains(dep{i}, "ft_getopt");  % Remove the dependency on ft_getopt, since ft_getopt is widely used
    dep{i} = dep{i}(indices & indices_dccnpath & indices_ft_getopt); 
  catch
    fprintf('cannot find "%s"\n', inlist{i});
    dep{i} = {};  % This function cannot be found, so no dependencies
  end
  num(i) = length(dep{i});
end

% Initialize a matrix to hold all dependencies
depmat = zeros(length(inlist), sum(num));

% Mark the dependencies with a value of 2
for i=1:length(inlist)
  offset = sum(num(1:(i-1)));
  depmat(i, (offset+1):(offset+num(i))) = 2;
end

% Drop functions with no dependencies
dep = dep(~num==0); 

% Define outlist and remove all double occurences
outlist=[dep{:}];
s = unique(outlist);
for i=1:length(s)
  sel = strmatch(s{i}, outlist);
  if ~isempty(sel)
      depmat(:,sel(1)) = max(depmat(:,sel),[], 2);
      outlist(sel(2:end)) = [];
      depmat(:,sel(2:end)) = [];
  end
end 