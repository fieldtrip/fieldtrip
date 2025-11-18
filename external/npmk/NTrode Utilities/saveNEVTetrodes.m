function saveNEVTetrodes(NEVFullFilename)

% saveNEVTetrodes
% 
% Opens saves a new NEV file and splits it into smaller NEV files
% containing only the tetrodes according to the associated CCF file. The
% CCF file should have the same name as the data file.
%
% NEVFullFilename:  The full path to the NEV to be opened.
%                   DEFAULT: If the filename is not provided, the user will
%                   be prompted to select a file.
%
% Use saveNEVTetrodes(NEVFullFilename)
%
%   Example: saveNEVTetrodes('c:\datafolder\datafile.nev');
%
%   The function will open datafile.nev and datafile.ccf and based on the
%   tetrode information saved in datafile.ccf, it will split datafile.nev
%   into smaller chunks that will only contain channels associated to that
%   particular tetrode.
%
%   Kian Torab
%   ktorab@blackrockmicro.com
%   Blackrock Microsystems
%   Version 1.0.0.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0:
%   - Initial release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Opening the file
if ~exist('NEVFullFilename', 'var')
    [dataFilename dataFolder] = getFile('*.nev');
    NEVFullFilename = [dataFolder dataFilename];
elseif exist(NEVFullFilename, 'file') ~= 2
        disp('The NEV file name does not exist.');
        return;
else
    [filepath filename fileext] = fileparts(NEVFullFilename);
    dataFolder = [filepath '/'];
    dataFilename = [filename fileext];
end
    

%% Openning the associated CCF file    
ccfFullFilename = [dataFolder dataFilename(1:end-3) 'ccf'];
if exist(ccfFullFilename, 'file') ~= 2 % for 2.x file type
    ccfFullFilename = [dataFolder dataFilename(1:end-8) '.ccf'];
    if exist(ccfFullFilename, 'file') ~= 2 % for TOC file type
        disp('Cannot find the associated CCF file.');
        disp('This function requires the CCF file used during the recording.');
        disp('The CCF must have the same name as the original recorded file.');
        return;
    end
end
ccf = openCCF(ccfFullFilename);

% Calculating the number of NTrode groups in the file
if isfield(ccf, 'NTrodeInfo')
    if isfield(ccf.NTrodeInfo, 'NTrodeMembers')
        numberOfNTrodeGroups = size(ccf.NTrodeInfo.NTrodeMembers,2);
    else
        disp('This data file does not contain any tetrodes.');
        return;
    end
else
    disp('There is an error in the CCF file.');
    return;
end

% Splitting the data file according to the information saved in the
% associated CCF file.
for idx = 1:numberOfNTrodeGroups
    fprintf('Saving tetrode %d containing channels %s.\n', idx, num2str(ccf.NTrodeInfo.NTrodeMembers{idx}));
    saveNEVSubSpikes(ccf.NTrodeInfo.NTrodeMembers{idx}, NEVFullFilename, 'tet');
end