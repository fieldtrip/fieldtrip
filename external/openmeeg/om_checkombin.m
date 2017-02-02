function om_checkombin
% Check if OpenMEEG binaries are installed and work.
%
% Copyright (C) 2010-2017, OpenMEEG developers

% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-30 11:15:51 +0200 (Thu, 30 Sep 2010) $
% $Revision$

[status,result] = system('om_assemble');
if status
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp('OpenMEEG binaries are not correctly installed')
    disp(' ')
    disp('Visit the OpenMEEG website for instructions')
    disp('http://openmeeg.github.io/')
    disp(' ')
    disp('See also the installation instructions on')
    disp('http://www.fieldtriptoolbox.org/faq/how_do_i_install_the_openmeeg_binaries')
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp(result);
    error('OpenMEEG not found')
end
