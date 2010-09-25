function filenames = read_dir(dirname, ext, prepend)
% READ_DIR Like the built-in dir command, but returns filenames as a cell array
% filenames = read_dir(dirname, ext)
%
% e.g., filenames = read_dir('images', '*.jpg')
% filenames{1} = 'foo.jpg', filenames{2} = 'foo2.jpg', etc
%
% read_dir(dirname, ext, 1) means pre-prend the directory name
%
% e.g., filenames = read_dir('images', '*.jpg', 1)
% filenames{1} = 'images/foo.jpg', filenames{2} = 'images/foo2.jpg', etc

if nargin < 3, prepend = 0; end

tmp = dir(fullfile(dirname, ext));
filenames = char(tmp(:).name); % filenames(f,:)
if prepend
  nfiles = size(filenames,1);
  if strncmp(computer,'PC',2)
    filenames = [repmat([dirname '\'], nfiles, 1) filenames];
  else
    filenames = [repmat([dirname '/'], nfiles, 1) filenames];
  end
end
filenames = num2cell(filenames, 2); % filenames{f}
