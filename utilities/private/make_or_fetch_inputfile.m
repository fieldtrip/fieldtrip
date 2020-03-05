function inputfile = make_or_fetch_inputfile(reproducescript_dir, filename, varname, value)

% MAKE_OR_FETCH_INPUTFILE is a helper function for ft_preamble_loadvar and ft_postamble_savevar, and
% is used for the cfg.reproducescript functionality.

% where to look for the hashes
hashfile = fullfile(reproducescript_dir, 'hashes.mat');
% use this filename in case nothing is found
inputfile = fullfile(reproducescript_dir, filename);

datahash = ft_hash(value);

if exist(hashfile, 'file')
  hashes = load(hashfile);
  fields = fieldnames(hashes);
  for k = 1:numel(fields)
    if strcmp(datahash, hashes.(fields{k}))
      % found a match, use the file given by this key
      inputfile = fields{k}(2:end); % strip off the leading 'f'
      inputfile = fullfile(reproducescript_dir, [inputfile '.mat']);
      return;
    end
  end
end

% if we're here, then there was no hash match; therefore we need to store
% the data into the inputfile
savevar(inputfile, varname, value, hashfile);

end