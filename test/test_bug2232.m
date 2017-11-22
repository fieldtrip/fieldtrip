function test_bug2232

% MEM 1gb
% WALLTIME 00:10:00

% TEST gifti
% TEST ft_read_headshape
% TEST ft_read_atlas

ft_hastoolbox('gifti', 1);

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2232');
d       = dir(datadir);
subdir  = {d.name}';
sel     = cellfun(@numel, subdir)>2;
subdir  = subdir(sel);

% check whether we can read *.surf.gii
for k = 1:numel(subdir)
  d = dir(fullfile(datadir,subdir{k},'*surf*'));
  if ~isempty(d)
    for m = 1:numel(d)
      fprintf('reading %s\n',d(m).name);
      bnd = ft_read_headshape(fullfile(datadir, subdir{k}, d(m).name));
      assert(isfield(bnd, 'pos') && isfield(bnd, 'tri'));
      
      hemi = strfind(d(m).name, '.L.');
      if ~isempty(hemi)
        hemi = '.L.';
      else
        hemi = '.R.';
      end
      dlabel = dir(fullfile(datadir,subdir{k},['*',hemi,'*label*']));
      if ~isempty(dlabel)
        for q = 1:numel(dlabel)
          fprintf('reading parcellation based on %s\n', dlabel(q).name);
          bnd = ft_read_atlas({fullfile(datadir,subdir{k},dlabel(q).name) fullfile(datadir,subdir{k},d(m).name)});
        end
      end
    end
  end
  
  % check whether there are subdirs
  d   = dir(fullfile(datadir,subdir{k}));
  sel = cellfun(@numel, {d.name}')>2;
  d   = d(sel);
  
  isdir     = find([d.isdir]);
  d         = d(isdir);
  if ~isempty(d)
    subsubdir = {d.name}';
    for p = 1:numel(subsubdir)
      d = dir(fullfile(datadir,subdir{k},subsubdir{p},'*surf*'));
      if ~isempty(d)
        for m = 1:numel(d)
          fprintf('reading %s\n',d(m).name);
          bnd = ft_read_headshape(fullfile(datadir, subdir{k}, subsubdir{p}, d(m).name));
          assert(isfield(bnd, 'pos') && isfield(bnd, 'tri'));
          
          hemi = strfind(d(m).name, '.L.');
          if ~isempty(hemi)
            hemi = '.L.';
          else
            hemi = '.R.';
          end
          dlabel = dir(fullfile(datadir,subdir{k},['*',hemi,'*label*']));
          if ~isempty(dlabel)
            for q = 1:numel(dlabel)
              fprintf('reading parcellation based on %s\n', dlabel(q).name);
              bnd = ft_read_atlas({fullfile(datadir,subdir{k},subsubdir{p},dlabel(q).name) fullfile(datadir,subdir{k},subsubdir{p},d(m).name)});
            end
          end

        end
      end
    end
  end
end

