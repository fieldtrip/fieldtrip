function codeKeeper = syncPatternDetectNEV(NEVStruct)

% syncPatternDetectNEV
%
% Finds the code for the SYNC pattern on the SYNC output of Blackrock
% Microsystems Neural Signal Processor (NSP) recorded by analog input 
% channel 16 and thresholded as extracted spikes and saved in the NEV file.
%
%   INPUT
%
%   NEVStruct:  The NEV structure containing the SYNC pulse.
%
%   OUTPUT
%
%   codeKeeper: A cell structure containing the unique code and their
%               respective timestamps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   USAGE EXAMPLE: 
%   
%   codeKeeper = syncPatternDetectNEV(NEVStruct)
%
%   In the example above, the NEV structure, containing the continueus 
%   signal from the SYNC pulse is passed to the function. The output will 
%   contain the unique codes parsed from the signal file and their 
%   respective timestamps.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   kian@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.0.0.0: July 7, 2014
%   - Initial release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
separationVar = 'D';

%% Detect the rising edges
edgeTS = NEVStruct.Data.Spikes.TimeStamp(NEVStruct.Data.Spikes.Electrode == 144);
differenceTS = double(diff(edgeTS));
pulseDifferenceLenght = mode(differenceTS);

convertedChar = [];
for idx = 1:length(differenceTS)
    if differenceTS(idx) < pulseDifferenceLenght*1.1;
        convertedChar = [convertedChar, '1'];
    elseif differenceTS(idx) < 9*pulseDifferenceLenght;
        convertedChar = [convertedChar, repmat('0', 1, round((differenceTS(idx) - pulseDifferenceLenght)/pulseDifferenceLenght)), '1'];
    else
        convertedChar = [convertedChar, separationVar];
    end
end

begTS = edgeTS(find(differenceTS > pulseDifferenceLenght * 12) + 1);
%% 
begTS = [edgeTS(1) begTS];

codeKeeper{1} = bin2dec(regexp(convertedChar, separationVar, 'split'))';
codeKeeper{2} = begTS;