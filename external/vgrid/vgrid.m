function executable = vgrid

% VGRID is a helper m-file to locate the binary executable

f = mfilename('fullpath');
[p, f, x] = fileparts(f);

% construct the full path to the binary executable
executable = fullfile(p, 'bin', lower(computer), 'vgrid');

if ispc
  executable = [executable '.exe'];
end

if ~exist(executable, 'file')
  error('the vgrid executable is not available for your platform, it was expected to be located at %s', executable);
end
