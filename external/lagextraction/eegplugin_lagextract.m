function eegplugin_extractlag( fig, try_strings, catch_strings)

% $Id: eegplugin_lagextract.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

uimenu( findobj(fig, 'tag', 'tools') , 'label', 'Extract lags', 'callback', ... 
           [try_strings.no_check '[EEG LASTCOM]= pop_extractlag(EEG);' catch_strings.store_and_hist '[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);' ]);
