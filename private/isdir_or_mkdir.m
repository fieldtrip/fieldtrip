function isdir_or_mkdir(p)

% ISDIR_OR_MKDIR Checks that a directory exists, or if not, creates the directory and
% all its parent directories.
%
% See also FOPEN_OR_ERROR

% add trailing file separator (slash or backslash)
p_full=[p filesep];

% we are looking for strings not containing the file separator, these are
% the subdirectories. We then store the last character position of each
% subdirectory in the variable 'match'.
pat=['[^\' filesep ']+'];
match=regexp(p_full, pat, 'end');

n_subdirs = numel(match);
for i=1:n_subdirs
  last_pos = match(i);
  d=p(1:last_pos);
  if ~isfolder(d)
    ft_notice('creating directory %s', d);
    mkdir(d);
  end
end
