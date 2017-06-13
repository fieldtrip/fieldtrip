function test_checkcode

% WALLTIME 00:20:00
% MEM 2gb

% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3309

%%

[v, p] = ft_version;
p = {p};
f = {};

while ~isempty(p)
  fprintf('looking in directory %s\n', p{1});
  nf = dir(p{1});
  p(1) = []; % remove this one
  for i=1:numel(nf)
    if nf(i).isdir && ~isequal(nf(i).name(1), '.') && ~isequal(nf(i).name, 'external')
      p{end+1} = fullfile(nf(i).folder, nf(i).name);
    else
      if ~nf(i).isdir && nf(i).name(end-1)=='.' && nf(i).name(end)=='m'
        f{end+1} = fullfile(nf(i).folder, nf(i).name);
      end
    end
  end
end

filelist = f;
clear f p

fprintf('found %d *.m files\n', numel(filelist));

%%

invalid = {
  'Use || instead of | as the OR operator in (scalar) conditional statements.'
  'Use && instead of & as the AND operator in (scalar) conditional statements.'
  };

m = {};
for i=1:numel(filelist)
  f = filelist{i};
  s = checkcode(f);
  m = union(m, {s.message});
  for j=1:numel(s)
    for k=1:numel(invalid)
      if strmatch(invalid{k},s(j).message)
        % pretty display
        fprintf('================================================================================\n');
        fprintf('ERROR: %s\n', invalid{k});
        fprintf('================================================================================\n');
        
        checkcode(f)
        error('please fix line %d in %s', s(j).line, f);
      end
    end
  end
end

% this shows the full list of warning messages
if false
  disp(m);
end

