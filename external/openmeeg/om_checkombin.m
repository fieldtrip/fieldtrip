function om_checkombin
% Check if OpenMEEG binaries are installed and work.
%
% Copyright (C) 2010, Alexandre Gramfort, INRIA

% $Id$
% $LastChangedBy$
% $LastChangedDate$
% $Revision$

[status,result] = system('om_assemble');
if status
    web('http://gforge.inria.fr/frs/?group_id=435')
    clc
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    disp('OpenMEEG binaries are not correctly installed')
    disp(' ')
    disp('Download OpenMEEG from')
    disp('http://gforge.inria.fr/frs/?group_id=435')
    disp(' ')
    disp('See wiki page for installation instructions:')
    disp('http://fieldtrip.fcdonders.nl/development/openmeeg/testinginstallation')
    disp('---------------------------------------------')
    disp('---------------------------------------------')
    error('OpenMEEG not found')
end
