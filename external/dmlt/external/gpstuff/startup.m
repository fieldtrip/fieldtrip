% If Matlab is started in GPstuff root directory, this startup.m
% adds subfolders automatically to the path
F = mfilename;
S = which(F);
if exist('OCTAVE_VERSION', 'builtin')
  subfolders={'diag' 'dist' 'gp' 'mc' 'misc' 'optim' 'xunit', 'octave_compat', 'inputparser'};
else
  subfolders={'diag' 'dist' 'gp' 'mc' 'misc' 'optim' 'xunit'};
end
for sf=subfolders
  addpath(strrep(S,[F '.m'],sf{:}))
end


% Alternatively copy following lines to startup.m in MATLAB startup folder
% and edit gpstuffroot to point where you have installed the GPstuff
%gpstuffroot='/...'
%addpath([gpstuffroot 'diag'])
%addpath([gpstuffroot 'dist'])
%addpath([gpstuffroot 'gp'])
%addpath([gpstuffroot 'mc'])
%addpath([gpstuffroot 'misc'])
%addpath([gpstuffroot 'optim'])
%addpath([gpstuffroot 'xunit'])

% If using Octave version of GPstuff, also add the following
%addpath([gpstuffroot 'inputparser'])
%addpath([gpstuffroot 'octave_compat'])