function failed_corrupt_matfiles

% MEM 8gb
% WALLTIME 10:00:00

files = findfiles(dccnpath('/home/common/matlab/fieldtrip/data/test'));
status = true(size(files));

for i=1:length(files)
  try
    tmp = load(files{i});
    clear tmp
    fprintf('passed reading %s\n', files{i});
  catch
    status(i) = false;
    warning('failed reading %s\n', files{i});
  end
end

if any(status==false)
  error('not all files could be read into MATLAB');
end

function files = findfiles(p)
d = dir(p);
d = d([d.bytes]<4*1e9); % skip files that are too large

files = {d(~[d.isdir]).name}';
files = files(~cellfun(@isempty, regexp(files, '\.mat$'))); % only select the mat files
files = files( cellfun(@isempty, regexp(files, '^\._')));   % exclude the OS X resource forks

for i=1:length(files)
  files{i} = fullfile(p, files{i});
end

dirs  = {d( [d.isdir]).name}';
dirs  = dirs(3:end); % the first two are . and ..

for i=1:length(dirs)
  files = cat(1, files, findfiles(fullfile(p, dirs{i})));
end

