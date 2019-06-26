function test_bug1067

% MEM 1500mb
% WALLTIME 00:10:00


% tests whether parameters that are in ft_freq* but not in ft_source* will
% be caught by ft_source* (and vice versa)

% change into fielctrip/test
cd(fileparts(mfilename('fullpath')));

% get all functions
source_functions = dir(fullfile('..','ft_source*.m'));
freq_functions = dir(fullfile('..','ft_freq*.m'));

% get function names
postfix = '.m';
source_prefix = 'ft_source';
source_fnc_names = cell(numel(source_functions), 1);
for i=1:numel(source_functions)
  source_fnc_names{i} = source_functions(i).name(numel(source_prefix)+1:end-numel(postfix));
end

freq_prefix = 'ft_freq';
freq_fnc_names = cell(numel(freq_functions), 1);
for i=1:numel(freq_functions)
  freq_fnc_names{i} = freq_functions(i).name(numel(freq_prefix)+1:end-numel(postfix));
end

% find common functions
freq_idx = find(ismember(freq_fnc_names, source_fnc_names));
source_idx = find(ismember(source_fnc_names, freq_fnc_names));
num_fncs = numel(freq_idx);

% loop through all functions
cfg_prefix = 'cfg.';
for i=1:num_fncs
  % open file content 
  s_content = cell2mat(textread(fullfile('..', source_functions(source_idx(i)).name), '%s', 'commentstyle', 'matlab', 'delimiter', '\n')');
  f_content = cell2mat(textread(fullfile('..', freq_functions(freq_idx(i)).name), '%s', 'commentstyle', 'matlab', 'delimiter', '\n')');
      
  % find cfg.* for each function
  s_cfg = strfind(s_content, cfg_prefix);
  f_cfg = strfind(f_content, cfg_prefix);
  
  % find options per function (hashmap like)
  s_opts = [];
  s_hash = [];
  delimiters = '=.%(,)\t '';/*-+{}<>[]~';
  for j=1:numel(s_cfg)
    s_opts = s_content(s_cfg(j)+numel(cfg_prefix):end);
    s_opts = textscan(s_opts, '%s', 'delimiter', delimiters, 'commentstyle', 'matlab');
    if isempty(s_opts{1}{1})
      continue;
    end
    s_hash.(s_opts{1}{1}) = true;
  end  

  f_opts = [];
  f_hash = [];
  for j=1:numel(f_cfg)
    f_opts = f_content(f_cfg(j)+numel(cfg_prefix):end);
    f_opts = textscan(f_opts, '%s', 'delimiter', delimiters, 'commentstyle', 'matlab');
    if isempty(f_opts{1}{1})
      continue;
    end
    f_hash.(f_opts{1}{1}) = true;
  end
  
  % find the not common options
  old_s_opts = fieldnames(s_hash);
  old_f_opts = fieldnames(f_hash);
  
  % options ft_source* uses but not ft_freq*
  s_opts = old_s_opts(~ismember(old_s_opts, old_f_opts));
  
  % options ft_freq* uses but not ft_source*
  f_opts = old_f_opts(~ismember(old_f_opts, old_s_opts));
  
  % output
  fprintf('\n');
  fprintf('Function %s uses the following fields that %s does not use:\n', ...
    source_functions(source_idx(i)).name, ...
    freq_functions(freq_idx(i)).name);
  for j=1:numel(s_opts)
    fprintf('\t cfg.%s\n', s_opts{j});
  end
 fprintf('\n');

  fprintf('Function %s uses the following fields that %s does not use:\n', ...
    freq_functions(freq_idx(i)).name, ...
    source_functions(source_idx(i)).name);
  for j=1:numel(f_opts)
    fprintf('\t cfg.%s\n', f_opts{j});
  end

end
