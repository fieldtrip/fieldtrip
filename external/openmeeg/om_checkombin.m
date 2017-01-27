function om_checkombin
% Check if OpenMEEG binaries are installed and work.
%
% Copyright (C) 2010, Alexandre Gramfort, INRIA

% Modification 2016 Nikolaas N. Oosterhof: to not use 'web' to open a web
% browser when OpenMEEG binaries are not correctly installed

% $Id$
% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-30 11:15:51 +0200 (Thu, 30 Sep 2010) $
% $Revision$

[status,result] = system('om_assemble');
if status
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp('OpenMEEG binaries are not correctly installed')
    disp(' ')
    disp('Download OpenMEEG from')
    disp('http://gforge.inria.fr/frs/?group_id=435')
    disp(' ')
    disp('See the OpenMEEG website at:');
    disp('http://openmeeg.gforge.inria.fr');
    disp(' ');
    disp('See the installation instructions on')
    disp('http://www.fieldtriptoolbox.org/faq/how_do_i_install_the_openmeeg_binaries')
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp(result);
    error('OpenMEEG not found')
end
