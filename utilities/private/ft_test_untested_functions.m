function untested_functions

% UNTESTED_FUNCTIONS finds FieldTrip high-level functions not tested by any test scripts
% 
% See also FIND_DEPENDENCY

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

% Get FieldTrip version and path
[ftver, ftpath] = ft_version;

% Define the folder containing test scripts
testFolder = fullfile(ftpath, 'test');

% List all files in the test folder
list = dir(testFolder);
testScripts  = {list.name};

% Exclude specific files and folders
excludeFiles = {'.', '..', 'test_eeglab_ft_integration_20140529T090217.mat', ...
'test_eeglab_ft_integration_20150528T150257.mat', ...
'test_issue1334.mat', ...
'test_spm_ft_integration.mat', ...
'trialfun_affcog.m', ...
'trialfun_stimon.m', ...
'trialfun_stimon_samples.m', 'invalid', 'private', 'README'};
validIndices = ~ismember({list.name}, excludeFiles);
testScripts  = testScripts(validIndices);

% Remove the ".m" extension from test script names
for i = 1:length(testScripts)
    testScripts{i} = strrep(testScripts{i}, '.m', '');
end

% Find the dependencies of all the test scripts
[dependencies, depmat] = find_dependency(testScripts);

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
disp('Untested functions:');
disp(char(untestedFunctions));
disp('Number of untested functions:');
disp(length(untestedFunctions));
