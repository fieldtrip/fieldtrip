function output = ntrodeGroups(ccf)

% ntrodeGroups
%
% This script takes in a CCF file and displays the information on the
% ntrode groups within the data file.
%
%   ccf:       Pass the CCF of interest. If no CCF is passed, the user will
%              be prompted to choose a CCF file.
%              DEFAULT: Will open Open File UI.
%    
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.1.0.0
%



if ~exist('ccf', 'var')
    ccf = openCCF;
end

for idx = 1:length(ccf.NTrodeInfo.NTrodeID)
    fprintf('NTrode group %d members: %s\n', idx, int2str(ccf.NTrodeInfo.NTrodeMembers{idx}));
end

output = ccf.NTrodeInfo.NTrodeMembers;
