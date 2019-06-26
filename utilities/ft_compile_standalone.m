function ft_compile_standalone(fieldtriproot)

% FT_COMPILE_STANDALONE compiles the FieldTrip functions along with
% the standalone entry function into a compiled executable.
%
% The compiled executable includes
%  - all main FieldTrip m-files
%  - all main FieldTrip m-files dependencies for as long as these
%    dependencies are in the fieldtrip modules and external toolboxes
%    on the path, MATLAB built-in, or toolbox/(stats/images/signal)
%    functions
%
% See also FT_STANDALONE, FT_COMPILE_MEX

% Copyright (C) 2011-2015 by Robert Oostenveld
%
% This file is part of FieldTrip.
%
% megconnectome is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% megconnectome is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with megconnectome.  If not, see <http://www.gnu.org/licenses/>.

% clear all variables, globals, functions and MEX links
% clear all; % don't do this because it clears the input arguments
clear global;
clear fun;
clear mex;

% v = ver('MATLAB');
% if ~strcmp(v.Version, '8.0')
%     error('the fieldtrip application should be compiled with MATLAB 2012b (8.0)');
% end

if nargin<1 || isempty(fieldtriproot)
  % this script is in fieldtrip/utilities
  fieldtriproot = fileparts(which(mfilename));
  fieldtriproot = fileparts(fieldtriproot);
end

origdir  = pwd;
builddir = fullfile(fieldtriproot, 'utilities');
bindir   = fullfile(fieldtriproot, 'bin');

cd(builddir);

% create a file with the timestamp of the compilation
fid = fopen('buildtimestamp.m', 'wt');
fprintf(fid, 'function s = buildtimestamp\n');
fprintf(fid, 's = ''%s'';\n', datestr(now));
fclose(fid);

if ~isfolder(bindir)
  mkdir(bindir);
end
cd(bindir);

fprintf('Using FieldTrip from "%s"\n', fieldtriproot);

% clean the path
restoredefaultpath;

%--------------------------------
% FIELDTRIP RELATED PATH SETTINGS

% add the path to fieldtrip
addpath(fieldtriproot);

% ensure that the path to the default modules is specified
clear ft_defaults;
ft_defaults;

% do not use my personal defaults, but rather FieldTrip standard defaults
global ft_default
ft_default = [];

% ensure that these special modules are also added
ft_hastoolbox('qsub', 1);
ft_hastoolbox('engine', 1);

% ensure that all desired external toolboxes are added to the path
% excluding spm2 and spm8 since they result in more than one spm on the path -> confusion in FT)
% excluding other toolboxes that are non-trivial to include in the compiled version
exclude = {
  '.'
  '..'
  '.svn'
  '.git'
  'eeglab'
  'dipoli'
  'dmlt'
  'iso2mesh'
  'simbio'
  'spm2'
  'spm8'
  'sqdproject'
  'yokogawa'
  'yokogawa_meg_reader'
  'ricoh_meg_reader'
  'egi_mff'     % java
  'mffmatlabio' % java
  'neuroscope'  % isstring
  };

extd = dir([fieldtriproot,'/external']);
extd = setdiff({extd([extd.isdir]).name}, exclude);
for k = 1:numel(extd)
  addpath(fullfile(fieldtriproot,'external', extd{k}));
end

%-------------------
% DO THE COMPILATION

cmd = ['mcc -R -singleCompThread -N -o fieldtrip -m ft_standalone ' ...
  ' -a ' fieldtriproot '/*.m' ...
  ' -a ' fieldtriproot '/utilities/*.m' ...
  ' -a ' fieldtriproot '/fileio/*.m' ...
  ' -a ' fieldtriproot '/forward/*.m' ...
  ' -a ' fieldtriproot '/inverse/*.m' ...
  ' -a ' fieldtriproot '/plotting/*.m' ...
  ' -a ' fieldtriproot '/statfun/*.m' ...
  ' -a ' fieldtriproot '/trialfun/*.m' ...
  ' -a ' fieldtriproot '/preproc/*.m' ...
  ' -a ' fieldtriproot '/qsub/*.m' ...
  ' -a ' fieldtriproot '/engine/*.m' ...
  ' -a ' fieldtriproot '/engine/*.m' ...
  ' -a ' fieldtriproot '/realtime/example/*.m' ...
  ' -a ' fieldtriproot '/realtime/online_eeg/*.m' ...
  ' -a ' fieldtriproot '/realtime/online_meg/*.m' ...
  ' -a ' fieldtriproot '/realtime/online_mri/*.m' ...
  ' -a ' fieldtriproot '/utilities/buildtimestamp.m' ...
  ' -p ' matlabroot    '/toolbox/signal' ...
  ' -p ' matlabroot    '/toolbox/images' ...
  ' -p ' matlabroot    '/toolbox/stats' ...
  ' -p ' matlabroot    '/toolbox/optim' ...
  ' -p ' matlabroot    '/toolbox/curvefit' ...
  ];
eval(cmd);

% somehow I don't manage to get this going with more than one directory to be added when calling mcc directly:
% mcc('-N', '-a', '/home/common/matlab/fieldtrip/*.m', '-o', fname, '-m', fname, '-p', [matlabroot,'/toolbox/signal:', matlabroot,'/toolbox/images:', matlabroot,'/toolbox/stats']);

% remove the additional files that were created during compilation
delete mccExcludedFiles.log
% delete readme.txt
% delete run_fieldtrip.sh
% delete buildtimestamp.m

fprintf('Finished compilation\n');
cd(origdir);
