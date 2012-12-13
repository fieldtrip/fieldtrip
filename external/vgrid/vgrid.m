function executable = vgrid

% VGRID is a helper m-file to locate the binary executable

f = mfilename('fullpath');
[p, f, x] = fileparts(f);

% construct the full path to the binary executable
executable = fullfile(p, 'bin', lower(computer), 'vgrid');

if ispc
  executable = [executable '.exe'];
end

