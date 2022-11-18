clear all
% this is the table that is currently read in from an existing file with data2bids
session_ids= {'ses-EcogLfpMedOn01';'ses-EcogLfpMedOn02'}
date=[(datetime('now'));char(datetime('now')+1)]
UPDRS = [ 7 ; 8]

sessions_tsv_1 = table(session_ids,date,UPDRS)

% this is what is given to the cfg of data2bids, here to update a UPDRS score
session_ids= {'ses-EcogLfpMedOn01'}
UPDRS=9
this = table(session_ids,UPDRS)

% the score cannot be updated because the datetime was not read in as a cell
mergetable(sessions_tsv_1, this, 'session_ids')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERROR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error using mergetable (line 159)
% The following error occurred converting from datetime to cell:
% Conversion to cell from datetime is not possible.
% 
% Error in data2bids (line 2178)
%       sessions_tsv = mergetable(sessions_tsv, this, 'session_id');
% 
% Error in BIDS_data_quality_control (line 104)
%     data2bids(cfg, intern_cfg.data);
%%%%%%%%%%%%%%%%%%%%%%%%%% SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% by changin the ft_read_tsv function to 'Text' instead of datetime, we get:
session_ids= {'ses-EcogLfpMedOn01';'ses-EcogLfpMedOn02'}
date={char(datetime('now'));char(datetime('now')+1)} %%% HERE IS THE DIFFERENCE
UPDRS = [ 7 ; 8]

sessions_tsv_1 = table(session_ids,date,UPDRS)

% same as above
session_ids= {'ses-EcogLfpMedOn01'}
UPDRS=9
this = table(session_ids,UPDRS)

% now it works and it writes out the sessions
mergetable(sessions_tsv_1, this, 'session_ids')

