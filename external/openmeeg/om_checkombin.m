function prefix = om_checkombin

% Check if OpenMEEG binaries are installed and work. This function
% will return the prefix required to execute the binaries.
%
% Copyright (C) 2010-2017, OpenMEEG developers

% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-30 11:15:51 +0200 (Thu, 30 Sep 2010) $
% $Revision$

if ispc
  [status, result] = system('om_assemble.exe');
  if status>0
    show_error(result);
  else
    prefix = '';
  end
  
else
  % the remainder of the code does not apply to windows
  
  location = {
    ''
    '/usr/bin'
    '/usr/local/bin'
    '/usr/local/openmeeg/bin'
    '/opt/bin'
    '/opt/openmeeg/bin'
    '/opt/openmeeg/2.4.1/bin'
    '/opt/openmeeg/2.4.0/bin'
    '/opt/openmeeg/2.3.0/bin'
    '/opt/openmeeg/2.2.0/bin'
    '/opt/openmeeg/2.1.0/bin'
    };
  
  % start with an empty return value
  prefix = '';
  
  for i=1:numel(location)
    % check whether the binary can be found
    [status, result] = system(sprintf('which %s', fullfile(location{i}, 'om_assemble')));
    if status==0
      prefix = location{i};
      % we found it
      break
    end
  end
  
  if status>0
    show_error(result);
  end
  
  % figure out where it is installed
  openmeeg_bin = fileparts(result);
  openmeeg_lib = fullfile(fileparts(openmeeg_bin), 'lib');
  
  % check whether the binary can be executed
  [status, result] = system([openmeeg_bin 'om_assemble']);
  
  % the failure might be due to the libraries not being found
  if status>0
    if ismac
      prefix = ['export DYLD_LIBRARY_PATH=' openmeeg_lib ' && ' openmeeg_bin filesep];
      [status, result] = system([prefix 'om_assemble']);
    elseif isunix
      prefix = ['export LD_LIBRARY_PATH=' openmeeg_lib ' && ' openmeeg_bin filesep];
      [status, result] = system([prefix 'om_assemble']);
    end
  end
  
  if status>0
    show_error(result);
  end
  
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
