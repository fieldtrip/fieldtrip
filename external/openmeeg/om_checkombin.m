function om_checkombin
% Check if OpenMEEG binaries are installed and work.
%
% Copyright (C) 2010, Alexandre Gramfort, INRIA

% $Id$
% $LastChangedBy: alegra $
% $LastChangedDate: 2010-09-30 11:15:51 +0200 (Thu, 30 Sep 2010) $
% $Revision$

[status,result] = system('om_assemble');
if status
    web('http://openmeeg.gforge.inria.fr')
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp('OpenMEEG binaries are not correctly installed')
    disp(' ')
    disp('Download OpenMEEG from')
    disp('http://gforge.inria.fr/frs/?group_id=435')
    disp(' ')
    disp('See the installation instructions on')
    disp('http://fieldtrip.fcdonders.nl/faq/how_do_i_install_the_openmeeg_binaries')
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp(result);
    error('OpenMEEG not found')
end
