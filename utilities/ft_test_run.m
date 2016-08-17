function status = ft_test_run(varargin)

% FT_TEST_RUN executes selected FieldTrip test scripts. It checks whether each test
% script runs without problems as indicated by an explicit error and posts the
% results on the FieldTrip dashboard.
%
% Use as
%   ft_test_run functionname
%
% Additional optional arguments are specified as key-value pairs and can include
%   dependency   = string
%   maxmem       = string
%   maxwalltime  = string
%
% Test functions should not require any input arguments.
% Output arguments of the test function will not be considered.
%
% See also FT_TEST_RESULT, FT_VERSION

% Copyright (C) 2016, Robert oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/donders/fieldtrip
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

optbeg = find(ismember(varargin, {'dependency', 'maxmem', 'maxwalltime'}));
if ~isempty(optbeg)
  optarg = varargin(optbeg:end);
  varargin = varargin(1:optbeg-1);
else
  optarg = {};
end

% get the optional input arguments
dependency  = ft_getopt(optarg, 'dependency', {});
maxmem      = ft_getopt(optarg, 'maxmem', inf);
maxwalltime = ft_getopt(optarg, 'maxwalltime', inf);

if ischar(dependency)
  % this should be a cell-array
  dependency = {dependency};
end

if ischar(maxwalltime)
  % it is probably formatted as HH:MM:SS
  maxwalltime = str2walltime(maxwalltime);
end

if ischar(maxmem)
  % it is probably formatted as XXmb, or XXgb, ...
  maxmem = str2mem(maxmem);
end

% get the version and the path
[revision, ftpath] = ft_version;

% testing a work-in-progress version is not supported
assert(istrue(ft_version('clean')), 'this requires all local changes to be committed');

%% determine the list of functions to test
if ~isempty(varargin) && exist(varargin{1}, 'file')
  functionlist = varargin;
else
  d = dir(fullfile(ftpath, 'test', 'test_*.m'));
  functionlist = {d.name}';
  for i=1:numel(functionlist)
    functionlist{i} = functionlist{i}(1:end-2); % remove the extension
  end
end

%% determine the list of files to test
filelist = cell(size(functionlist));
for i=1:numel(functionlist)
  filelist{i} = which(functionlist{i});
end

fprintf('considering %d test scripts for execution\n', numel(filelist));

%% make a subselection based on the filters
sel = true(size(filelist));
mem = zeros(size(filelist));
tim = zeros(size(filelist));

for i=1:numel(filelist)
  fid = fopen(filelist{i}, 'rt');
  str = fread(fid, [1 inf], 'char=>char');
  fclose(fid);
  line = tokenize(str, 10);
  
  if ~isempty(dependency)
    sel(i) = false;
  else
    sel(i) = true;
  end
  
  for k=1:numel(line)
    for j=1:numel(dependency)
      [s, e] = regexp(line{k}, sprintf('%% TEST.*%s.*', dependency{j}), 'once', 'start', 'end');
      if ~isempty(s)
        sel(i) = true;
      end
    end
    
    [s, e] = regexp(line{k}, '% WALLTIME.*', 'once', 'start', 'end');
    if ~isempty(s)
      s = s + length('% WALLTIME'); % strip this part
      tim(i) = str2walltime(line{k}(s:e));
    end
    
    [s, e] = regexp(line{k}, '% MEM.*', 'once', 'start', 'end');
    if ~isempty(s)
      s = s + length('% MEM'); % strip this part
      mem(i) = str2mem(line{k}(s:e));
    end
  end % for each line
  
end % for each function/file

fprintf('%3d scripts do not meet the requirements for dependencies\n', sum(~sel));
fprintf('%3d scripts do not meet the requirements for memory\n',       sum(mem>maxmem));
fprintf('%3d scripts do not meet the requirements for walltime \n',    sum(tim>maxwalltime));

% remove test scripts that exceed walltime or memory
sel(tim>maxwalltime) = false;
sel(mem>maxmem)      = false;

% make the subselection of functions to test
functionlist = functionlist(sel);

fprintf('executing %d test scripts\n', numel(functionlist));

%% run over all tests
for i=1:numel(functionlist)
  
  close all
  fprintf('================================================================================\n');;
  fprintf('=== evaluating %s\n', functionlist{i});
  try
    stopwatch = tic;
    eval(functionlist{i});
    status = true;
    runtime = round(toc(stopwatch));
    fprintf('=== %s PASSED in %d seconds\n', functionlist{i}, runtime);
  catch
    status = false;
    runtime = round(toc(stopwatch));
    fprintf('=== %s FAILED in %d seconds\n', functionlist{i}, runtime);
  end
  close all
  
  result = [];
  result.matlabversion    = version('-release');
  result.fieldtripversion = revision;
  result.branch           = ft_version('branch');
  result.hostname         = gethostname;
  result.user             = getusername;
  result.result           = status;
  result.runtime          = runtime;
  result.functionname     = functionlist{i};
  
  options = weboptions('MediaType','application/json');
  webwrite('http://dashboard.fieldtriptoolbox.org/api', result, options);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function walltime = str2walltime(str)
str = lower(strtrim(str));
walltime = str2double(str);
if isnan(walltime)
  str = strtrim(str);
  hms = sscanf(str, '%d:%d:%d'); % hours, minutes, seconds
  walltime = 60*60*hms(1) + 60*hms(2) + hms(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mem = str2mem(str)
str = lower(strtrim(str));
mem = str2double(str);
if isnan(mem)
  if ~isempty(regexp(str, '[0-9]*kb', 'once'))
    mem = str2double(str(1:end-2)) * 2^10;
  elseif ~isempty(regexp(str, '[0-9]*mb', 'once'))
    mem = str2double(str(1:end-2)) * 2^20;
  elseif ~isempty(regexp(str, '[0-9]*gb', 'once'))
    mem = str2double(str(1:end-2)) * 2^30;
  elseif ~isempty(regexp(str, '[0-9]*tb', 'once'))
    mem = str2double(str(1:end-2)) * 2^40;
  else
    mem = str2double(str);
  end
end
