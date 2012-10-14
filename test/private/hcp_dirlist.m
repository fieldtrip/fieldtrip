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
