function ft_test_compare(varargin)

% FT_TEST_COMPARE

command  = varargin{1};
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

% construct the query string
query = '?';
queryparam = {'matlabversion', 'fieldtripversion', 'hostname', 'user', 'branch'};
for i=1:numel(queryparam)
  val = ft_getopt(optarg, queryparam{i});
  if ~isempty(val)
    query = [query sprintf('%s=%s&', queryparam{i}, val)];
  end
end

options = weboptions('ContentType','json'); % this returns the results as MATLAB structure

arg1 = varargin{1};
arg2 = varargin{2};

switch command
  case 'comparerevision'
    functionname = {};
    dashboard    = {};
    for i=1:numel(varargin)
      dashboard{i} = webread(['http://dashboard.fieldtriptoolbox.org/api/' query sprintf('&fieldtripversion=%s', varargin{i})], options);
      assert(~isempty(dashboard{i}), 'no tests were returned for revision %d', i);
      extraname    = {dashboard{i}.functionname};
      functionname = cat(1, functionname(:), extraname(:));
    end
    functionname = unique(functionname);
    
    % represent a summary of the results in a struct-array
    results = struct();
    for i=1:numel(functionname)
      results(i).function = functionname{i};
      for j=1:numel(varargin)
        sel = find(strcmp({dashboard{j}.functionname}, functionname{i}));
        results(i).(['rev_' varargin{j}]) = getresult(dashboard{j}, sel);
      end % for each functionname
    end % for each revision
    
    % convert the struct-array to a table
    table = struct2table(results);
    fprintf('%s\n', table{:});
    
  otherwise
    error('unsupported command "%s"', command);
end % switch command

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = padto(str, n)
if n>length(str)
  str = [str repmat(' ', [1 n-length(str)])];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str = getresult(dashboard, sel)
if isempty(sel)
  str = 'missing';
elseif all(istrue([dashboard(sel).result]))
  str = 'passed';
elseif all(~istrue([dashboard(sel).result]))
  str = 'failed';
else
  str = 'ambiguous';
end

