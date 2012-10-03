function test_ft_datatype

% TEST test_ft_datatype
% TEST ft_datatype ft_datatype_comp ft_datatype_mvar ft_datatype_source ft_datatype_dip ft_datatype_parcellation ft_datatype_spike ft_datatype_freq ft_datatype_raw ft_datatype_timelock ft_datatype_headmodel ft_datatype_segmentation ft_datatype_volume ft_datatype ft_datatype_sens

dirlist = {
  '/home/common/matlab/fieldtrip/data/test/latest'
  '/home/common/matlab/fieldtrip/data/test/20111231'
  '/home/common/matlab/fieldtrip/data/test/20110630'
  '/home/common/matlab/fieldtrip/data/test/20101231'
  '/home/common/matlab/fieldtrip/data/test/20100630'
  '/home/common/matlab/fieldtrip/data/test/20091231'
  '/home/common/matlab/fieldtrip/data/test/20090630'
  '/home/common/matlab/fieldtrip/data/test/20081231'
  '/home/common/matlab/fieldtrip/data/test/20080630'
  '/home/common/matlab/fieldtrip/data/test/20071231'
  '/home/common/matlab/fieldtrip/data/test/20070630'
  '/home/common/matlab/fieldtrip/data/test/20061231'
  '/home/common/matlab/fieldtrip/data/test/20060630'
  '/home/common/matlab/fieldtrip/data/test/20051231'
  '/home/common/matlab/fieldtrip/data/test/20050630'
  '/home/common/matlab/fieldtrip/data/test/20040623'
  '/home/common/matlab/fieldtrip/data/test/20031128'
  };

for j=1:length(dirlist)
  filelist = hcp_filelist(dirlist{j});
  
  [~, ~, x] = cellfun(@fileparts, filelist, 'uniformoutput', false);
  sel = strcmp(x, '.mat');
  filelist = filelist(sel);
  clear p f x
  
  for i=1:length(filelist)
    
    try
      fprintf('processing data structure %d from %d\n', i, length(filelist));
      var = loadvar(filelist{i});
      disp(var)
    catch
      % some of the mat files are corrupt, this should not spoil the test
      disp(lasterr);
      continue
    end
    
    type = 'unknown';
    
    if ~isempty(regexp(filelist{i}, '/raw/'))
      type = 'raw';
    elseif ~isempty(regexp(filelist{i}, '/comp/'))
      type = 'comp';
    elseif ~isempty(regexp(filelist{i}, '/timelock/'))
      type = 'timelock';
    elseif ~isempty(regexp(filelist{i}, '/freq/'))
      type = 'freq';
    elseif ~isempty(regexp(filelist{i}, '/source/'))
      type = 'source';
    elseif ~isempty(regexp(filelist{i}, '/volume/')) || ~isempty(regexp(filelist{i}, '/mri/'))
      type = 'volume';
    end
    
    switch type
      case 'raw'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'comp'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
        type = 'raw'; % comp data is a special type of raw data
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'timelock'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'freq'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'source'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      case 'volume'
        assert(ft_datatype(var, type), sprintf('%s did not contain %s data', filelist{i}, type));
      otherwise
        warning('not testing %s', filelist{i});
        % do nothing
    end % switch
  end % for filelist
end % for dirlist


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [list, numdirs, numfiles] = hcp_dirlist(basedir, recursive)

if nargin<2
  recursive = true;
end

if ~isdir(basedir)
  error('directory "%s" does not exist', basedir)
end

list = dir(basedir);

% remove all non-directories and hidden directories
list = list([list.isdir]);
hidden = false(size(list));
for i=1:length(list)
  hidden(i) = list(i).name(1)=='.';
end
list = list(~hidden);

% convert to cell-array
list = {list.name};
list = list(:);
for i=1:length(list)
  list{i} = fullfile(basedir, list{i});
end

list = sort(list);
numdirs = nan(size(list));
numfiles = nan(size(list));

for i=1:length(list)
  content = dir(list{i});
  numdirs(i) = sum([content.isdir]) - 2;
  numfiles(i) = length(content) - numdirs(i) - 2;
end

if recursive
  sub_list = cell(size(list));
  sub_numdirs = cell(size(list));
  sub_numfiles = cell(size(list));
  for i=1:length(list)
    [sub_list{i}, sub_numdirs{i}, sub_numfiles{i}] = hcp_dirlist(list{i}, recursive);
  end
  list = cat(1, list, sub_list{:});
  numdirs = cat(1, numdirs, sub_numdirs{:});
  numfiles = cat(1, numfiles, sub_numfiles{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function list = hcp_filelist(basedir)

dirlist = hcp_dirlist(basedir, true);
dirlist{end+1} = basedir;
list = {};

for i=1:length(dirlist)
  f = dir(dirlist{i});
  f = f(~[f.isdir]);
  f = {f.name};
  for j=1:length(f)
    f{j} = fullfile(dirlist{i}, f{j});
  end
  list = cat(1, list, f(:));
end

list = sort(list);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = loadvar(filename, varname)

if nargin<2
  fprintf('reading variable from file ''%s''\n', filename);
else
  fprintf('reading ''%s'' from file ''%s''\n', varname, filename);
end

var = whos('-file', filename);

if length(var)==1
  filecontent = load(filename); % read the one variable in the file, regardless of how it is called
  value = filecontent.(var.name);
  clear filecontent
else
  filecontent = load(filename, varname);
  value = filecontent.(varname); % read the variable named according to the input specification
  clear filecontent
end

