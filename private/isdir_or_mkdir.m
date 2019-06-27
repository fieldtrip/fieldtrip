function isdir_or_mkdir(p)

% ISDIR_OR_MKDIR Checks that a directory exists, or if not, creates the directory and
% all its parent directories.
%
% See also FOPEN_OR_ERROR

tok = tokenize(p, filesep);
if isempty(tok{1})
  tok{1} = filesep;
end
for i=1:numel(tok)
  d = fullfile(tok{1:i});
  if ~isfolder(d)
    ft_notice('creating directory %s', d);
    mkdir(d);
  end
end
