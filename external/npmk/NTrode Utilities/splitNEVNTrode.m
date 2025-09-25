function splitNEVNTrode

% splitNEVNTrode
% 
% Opens and splits an NEV file based on its NTrode groups. It depends on
% openCCF file.
%
% Use splitNEVNTrode
% 
% This function does not take any inputs.
%
%   Example 1: 
%   splitNEVNTrode;
%
%   In the example above, the user will be prompted to select a CCF file
%   first. The CCF contains the ntrode grouping infromation. Then the user
%   will be prompted to select a NEV file. The script will then split the
%   NEV file into smaller NEV files containing channels in given ntrode
%   groups. For example, if ntrode group one consists of channels 1,3,5,
%   and 12, then using splitNEVNTrode will split the file into a smaller
%   NEV files that contains those channels only. If there are multiple
%   ntrodes then the files will split into multiple smaller files, equal in
%   number of the ntrode groups.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: January 18, 2016
%   - Initial release.
%
% 1.1.0.0: October 10, 2016 - Saman Hagh Gooie
%   - Bug fixes with file loading 
%   - Fixed the file extension used for saving 
%
% 1.2.0.0: October 27, 2020
%   - Removed junk characters from the file.
%
% 1.2.1.0: November 17, 2020
%   - Minor bug fixes and general code clean up - @David Kluger
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validating input parameter
ccf = openCCF;
splitCount = unique(ccf.NTrodeInfo.NTrodeID);
splitCount = splitCount(splitCount>0);

% Getting the file name

if ~ismac
    [fname, path] = getFile('*.nev', 'Choose an NEV file...');
else
    [fname, path] = getFile('*.nev', 'Choose an NEV file...');
end

if fname == 0
    disp('No file was selected.');
    if nargout
        clear variables;
    end
    return;
end
    
% Loading the file    
for idx = splitCount
    % Determining whether tetrode channel is recorded and valid in NSx
    tetrodeChannels = ccf.NTrodeInfo.NTrodeMembers{idx};
    NEV = openNEV([path, fname(1:end-4) '.nev'], ['c:' num2str(tetrodeChannels)]); % modified by SH 05-Oct-2016
    newFileName = [path fname(1:end-4) '-tet' sprintf('%03d', idx) fname(end-3:end)];
    saveNEV(NEV, newFileName, 'noreport');
end