function [result] = ft_test_run(varargin)

% FT_TEST_RUN

% Copyright (C) 2017, Robert Oostenveld
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

narginchk(1, inf);
command = varargin{1};
assert(isequal(command, 'run') || isequal(command, 'inventorize'));
varargin = varargin(2:end);

optbeg = find(ismember(varargin, {'dependency', 'dccnpath', 'maxmem', 'maxwalltime', 'upload', 'sort', 'returnerror'}));
if ~isempty(optbeg)
  optarg   = varargin(optbeg:end);
  varargin = varargin(1:optbeg-1);
else
  optarg = {};
end

% varargin contains the file (or files) to test
% optarg contains the command-specific options

% get the optional input arguments
dependency  = ft_getopt(optarg, 'dependency', {});
hasdccnpath = ft_getopt(optarg, 'dccnpath');  % default is handled below
maxmem      = ft_getopt(optarg, 'maxmem', inf);
maxwalltime = ft_getopt(optarg, 'maxwalltime', inf);
upload      = ft_getopt(optarg, 'upload', 'yes'); % this will be set to 'no' in case FieldTrip version is not clean
sortarg     = ft_getopt(optarg, 'sort', 'alphabetical');
returnerror = ft_getopt(optarg, 'returnerror', 'no');

if ischar(dependency)
  % this should be a cell-array
  dependency = {dependency};
end

if isempty(hasdccnpath)
  % true when central storage is available, false otherwise
  hasdccnpath = exist(dccnpath('/home/common/matlab/fieldtrip/data/'), 'dir');
end

if ischar(maxwalltime)
  % it is probably formatted as HH:MM:SS, convert to seconds
  maxwalltime = str2walltime(maxwalltime);
end

if ischar(maxmem)
  % it is probably formatted as XXmb, or XXgb, convert to bytes
  maxmem = str2mem(maxmem);
end

% get the version and the path
[revision, ftpath] = ft_version;

%% determine the list of functions to test
if ~isempty(varargin) && exist(varargin{1}, 'file')
  functionlist = varargin;
else
  addpath(fullfile(ftpath, 'test'));
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

fprintf('considering %d test scripts\n', numel(filelist));

%% make a subselection based on the filters
dep  = true(size(filelist));
file = false(size(filelist)); % by default ignore file reads
mem  = zeros(size(filelist));
tim  = zeros(size(filelist));

for i=1:numel(filelist)
  fid = fopen(filelist{i}, 'rt');
  str = fread(fid, [1 inf], 'char=>char');
  fclose(fid);
  line = tokenize(str, 10);
  
  if ~isempty(dependency)
    dep(i) = false;
  else
    dep(i) = true;
  end
  
  for k=1:numel(line)
    for j=1:numel(dependency)
      % search for the dependencies in each of the test functions
      [s, e] = regexp(line{k}, sprintf('%% DEPENDENCY.*%s.*', dependency{j}), 'once', 'start', 'end');
      if ~isempty(s)
        dep(i) = true;
      end
    end
    
    if ~istrue(hasdccnpath)
      % search for the occurence of the DCCNPATH function in each of the test functions
      [s, e] = regexp(line{k}, 'dccnpath', 'once', 'start', 'end');
      if ~isempty(s)
        file(i) = true;
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

fprintf('%3d scripts are excluded due to the dependencies\n',                  sum(~dep));
fprintf('%3d scripts are excluded due to loading files from the DCCN path\n',  sum(file));
fprintf('%3d scripts are excluded due to the requirements for memory\n',       sum(mem>maxmem));
fprintf('%3d scripts are excluded due to the requirements for walltime \n',    sum(tim>maxwalltime));

% make selection of test scripts, remove scripts that exceed walltime or memory
sel                  = true(size(filelist));
sel(tim>maxwalltime) = false;
sel(mem>maxmem)      = false;
sel(dep==false)      = false;
sel(file==true)      = false;

% make the subselection of functions to test
functionlist = functionlist(sel);
tim = tim(sel);
mem = mem(sel);

switch (sortarg)
  case {'alphabet' 'alphabetic' 'alphabetical'}
    [functionlist, indx] = sort(functionlist);
    tim = tim(indx);
    mem = mem(indx);
  case {'mem' 'memory'}
    [mem, indx]  = sort(mem);
    tim          = tim(indx);
    functionlist = functionlist(indx);
  case {'walltime' 'time'}
    [tim, indx]  = sort(tim);
    mem          = mem(indx);
    functionlist = functionlist(indx);
  case 'random'
    indx = randperm(numel(functionlist));
    functionlist = functionlist(indx);
    tim          = tim(indx);
    mem          = mem(indx);
  otherwise
    ft_error('incorrect specification for ''sort''');
end

% this will hold all results
result = [];

if strcmp(command, 'inventorize')
  % do not execute anything, only return the function names
  for i=1:numel(functionlist)
    result(i).functionname = functionlist{i};
    result(i).mem          = mem(i);
    result(i).walltime     = tim(i);
  end
  fprintf('selecting %d test scripts\n', numel(functionlist));
  return
else
  fprintf('executing %d test scripts\n', numel(functionlist));
end

% uploading results to the dashboard from a work-in-progress version is not supported
if istrue(upload) && ~istrue(ft_version('clean'))
  warning('FieldTrip version is not clean, not uploading results to the dashboard')
  upload = 'no';
end

%% run over all tests
for i=1:numel(functionlist)
  
  close all
  fprintf('================================================================================\n');;
  fprintf('=== evaluating %s\n', functionlist{i});
  try
    stopwatch = tic;
    eval(functionlist{i});
    passed = true;
    runtime = round(toc(stopwatch));
    fprintf('=== %s PASSED in %d seconds\n', functionlist{i}, runtime);
  catch me
    passed = false;
    runtime = round(toc(stopwatch));
    % show the error with the stack trace
    fprintf('Error using %s (line %d)\n', me.stack(1).name, me.stack(1).line);
    disp(me.message)
    for j=2:numel(me.stack)
      fprintf('Error in %s (line %d)\n', me.stack(j).name, me.stack(j).line);
    end
    fprintf('=== %s FAILED in %d seconds\n', functionlist{i}, runtime);
  end
  close all
  
  result(i).matlabversion    = version('-release');
  result(i).fieldtripversion = revision;
  result(i).branch           = ft_version('branch');
  result(i).arch             = computer('arch');
  result(i).hostname         = gethostname;
  result(i).user             = getusername;
  result(i).passed           = passed;
  result(i).runtime          = runtime;
  result(i).functionname     = functionlist{i};
  result(i).mem              = mem(i);
  result(i).walltime         = tim(i);
  
  if istrue(upload)
    try
      % the weboptions function is available in 2014b onward, but behaves inconsistently
      options = weboptions('MediaType','application/json');
    catch
      options = [];
    end
    url = 'http://dashboard.fieldtriptoolbox.org/api/';
    webwrite(url, result(i), options);
  end
  
  if strcmp(returnerror, 'immediate') && ~passed
    % report only on the tests completed sofar and give an error
    printstruct_as_table(result);
    rethrow(me);
  end
  
end

if strcmp(returnerror, 'final') && any(~[result.passed])
  % report on the complete batch and give an error
  printstruct_as_table(result);
  failed = find(~[result.passed]);
  error('%d of the test scripts failed', numel(failed));
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
