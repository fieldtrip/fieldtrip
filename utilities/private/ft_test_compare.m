function ft_test_compare(varargin)

% FT_TEST_COMPARE

assert(numel(varargin)>1, 'not enough input arguments')
command = varargin{1};
feature = varargin{2};
assert(isequal(command, 'compare'));
varargin = varargin(3:end);

optbeg = find(ismember(varargin, {'matlabversion', 'fieldtripversion', 'user', 'hostname', 'branch'}));
if ~isempty(optbeg)
  optarg   = varargin(optbeg:end);
  varargin = varargin(1:optbeg-1);
else
  optarg = {};
end

% varargin contains the file (or files) to test
% optarg contains the command-specific options

% construct the query string that will be passed in the URL
query = '?';
queryparam = {'matlabversion', 'fieldtripversion', 'hostname', 'user', 'branch'};
for i=1:numel(queryparam)
  val = ft_getopt(optarg, queryparam{i});
  if ~isempty(val)
    query = [query sprintf('%s=%s&', queryparam{i}, val)];
  end
end

options = weboptions('ContentType','json'); % this returns the results as MATLAB structure

results = cell(size(varargin));
functionname = {};
for i=1:numel(varargin)
  results{i} = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&%s=%s', feature, varargin{i})], options);
  assert(~isempty(results{i}), 'no results were returned for %s %s', feature, varargin{i});
  functionname = cat(1, functionname(:), {results{i}.functionname}');
end

% find the joint set of all functions
functionname = unique(functionname);

% represent a summary of all results in a struct-array
summary = struct();
for i=1:numel(functionname)
  summary(i).function = functionname{i};
  for j=1:numel(varargin)
    sel = find(strcmp({results{j}.functionname}, functionname{i}));
    summary(i).(fixname(varargin{j})) = getresult(results{j}, sel);
  end % for each functionname
end % for each of the features

% convert the struct-array to a table
table = struct2table(summary);
fprintf('%s\n', table{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = fixname(str)
str = ['x_' str];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = getresult(result, index)
if isempty(index)
  str = 'missing';
elseif all(istrue([result(index).result]))
  str = 'passed';
elseif all(~istrue([result(index).result]))
  str = 'failed';
else
  str = 'ambiguous';
end

