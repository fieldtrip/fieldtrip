function prefix = om_checkombin

% Check if OpenMEEG binaries are installed and work.
%
% Copyright (C) 2010-2017, OpenMEEG developers

% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-30 11:15:51 +0200 (Thu, 30 Sep 2010) $
% $Revision$

% start with an empty return value
prefix = '';

% check whether the binary is on the path
[status, result] = system('which om_assemble');
if status>0
  show_error(result);
end

% figure out where it is installed
openmeeg_bin = fileparts(result);
openmeeg_lib = fullfile(fileparts(openmeeg_bin), 'lib');

% check whether the binary can be executed
[status, result] = system('om_assemble');

% the failure might be due to the libraries not being found
if status>0
  if ismac
    prefix = ['export DYLD_LIBRARY_PATH=' openmeeg_lib ' && '];
    [status, result] = system([prefix 'om_assemble']);
  elseif isunix
    prefix = ['export LD_LIBRARY_PATH=' openmeeg_lib ' && '];
    [status, result] = system([prefix 'om_assemble']);
  elseif ispc
    % don't know how to determine this on Windows
  end
end

if status>0
  show_error(result);
end

if nargout==0
  clear prefix
end


function show_error(result)
disp('==============================================================================');
disp('OpenMEEG binaries are not correctly installed')
disp(' ')
disp('Visit the OpenMEEG website for instructions')
disp('http://openmeeg.github.io/')
disp(' ')
disp('See also the installation instructions on')
disp('http://www.fieldtriptoolbox.org/faq/how_do_i_install_the_openmeeg_binaries')
disp('==============================================================================');
disp(result);
error('OpenMEEG not found')
