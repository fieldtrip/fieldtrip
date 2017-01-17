function results = ft_test_report(varargin)

% FT_TEST_REPORT

assert(numel(varargin)>0, 'not enough input arguments')
command  = varargin{1};
assert(isequal(command, 'report'));
varargin = varargin(2:end);

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

results = webread(['http://dashboard.fieldtriptoolbox.org/api/' query], options);

% remove some of the fields
results = removefields(results, {'x_id', 'createDate'});

% convert the struct-array to a table
table = struct2table(results);
fprintf('%s\n', table{:});
