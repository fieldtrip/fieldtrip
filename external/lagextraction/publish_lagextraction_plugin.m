% $Id: publish_lagextraction_plugin.m 4 2009-08-15 21:10:35Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-08-15 17:10:35 -0400 (Sam, 15 aoû 2009) $
% $Revision: 4 $

clear
opts.stylesheet = 'doc/agramfort.xsl';
opts.outputDir = 'doc';
opts.maxWidth = 600;
opts.catchError = false;
opts.evalCode = true;
% opts.evalCode = false;

% publish('lagextraction_demo',opts);
publish('lagextraction_tutorial',opts);
zip('../lagextraction_eeglab_plugin.zip','../lagextraction');
