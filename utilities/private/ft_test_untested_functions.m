function ft_test_untested_functions(varargin)

% FT_TEST_UNTESTED_FUNCTIONS documentation is included inside ft_test
% documentation.
% 
% See also FT_TEST

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

narginchk(1, 1); % Check that there is exactly one input argument
command = varargin{1};
assert(isequal(command, 'untested_functions'));

% Get FieldTrip version and path
[ftver, ftpath] = ft_version;

% Define the folder containing test scripts
testFolder = fullfile(ftpath, 'test');

% List all files in the test folder that start with test_* or inspect_*
list1 = dir(fullfile(testFolder,'test_*'));
list2 = dir(fullfile(testFolder,'inspect_*'));
list  = [list1; list2];
testScripts  = {list.name};

% Remove the ".m" extension from test script names
for i = 1:length(testScripts)
    testScripts{i} = strrep(testScripts{i}, '.m', '');
end

% Find the dependencies of all the test scripts
[dependencies, depmat] = ft_test_find_dependency('untested_functions', testScripts(1:10));

list    = dir(ftpath);
indices = endsWith({list.name}, {'.m'});
list    = list(indices); % List all .m files in the FieldTrip path 
indices = ~contains({list.name},  {'Contents.m'});
list    = list(indices); % Excluding "Contents.m"

% Remove the ".m" extension from function names
functionFiles = {list.name};
for i = 1:length(functionFiles)
    functionFiles{i} = strrep(functionFiles{i}, '.m', '');
end

% Replace dependencies with the filename, excluding the path
for i = 1:length(dependencies)
    [p, f, x] = fileparts(dependencies{i});
    dependencies{i} = f;
end

% Compare functionFiles with dependencies and find what's missing
testedFunctionsIndices = zeros(1, length(functionFiles));
for i = 1:length(dependencies)
    isDependency = contains(functionFiles, dependencies{i});
    testedFunctionsIndices = isDependency | testedFunctionsIndices;
end

testedFunctions   = functionFiles(testedFunctionsIndices);
untestedFunctions = setdiff(functionFiles, testedFunctions);

% Display the untested functions
fprintf('\n ------------------- \n')
disp('Untested functions:');
disp(char(untestedFunctions));
fprintf('\n Number of untested functions: %d \n\n', length(untestedFunctions))

