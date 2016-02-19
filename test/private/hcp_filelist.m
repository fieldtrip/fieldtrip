
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
