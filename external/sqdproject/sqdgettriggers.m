function [triggerLinesGood, triggerCounts, triggersList, fs] = sqdgettriggers(sqdfile, triggerLinesPotential, triggerHighThresh, triggerLowThresh)

%SQDGETTRIGGERS Detects trigger lines and triggers on those lines from an MEG sqd file.
%
%    [triggerLines, triggerCounts, triggersList, fs] = SQDGETTRIGGERS(sqdfile, triggerLinesPotential, triggerHighThresh, triggerLowThresh);
%
%    OUTPUT variables:
%    triggerLines: which lines were actually used as trigger lines
%    triggerCounts: how many triggers were detected on each trigger line
%    triggersList: a cell array (each cell corresponds to a
%       different trigger line) in which each cell is a list of the trigger
%       events. The trigger events are labeled by sample number not time. To
%       convert into time, divide the sample numbers by the sample
%       frequency.
%    fs: the sample frequency, in Hz.
%
%    INPUT variables:
%    sqdfile: the MEG file holding the trigger (and MEG) data
%    triggerLinesPotential: a list of which lines should be inspected for
%       triggers. The default is [160:191], i.e. all channels not containing
%       magnetic field data. This function is much more efficient, though, if
%       the triggerLinesPotential is the list of channels actually used as
%       trigger lines. Because of the way the sqd file is read, however, the
%       fastest method is to set triggerLinesPotential to be list of all
%       channels from the lowest used to the highest used.
%    triggerHighThresh (rarely needed): the threshold, in Volts, above
%       which a trigger line goes high. This value is nominally slightly
%       larger than 5 V, so the default value is 5 V.
%    triggerHighThresh (rarely needed): the threshold, in Volts, below
%       which a trigger line is considered to be low. This
%       value is nominally 0 V, so the default value is 0.2 V.

%
%    Version 0.9 beta 2
%    9 January 2010 
%    
%    by Jonathan Z. Simon


% Create triggerLinesPotential if not supplied or if equal to []. The
%    The default is all high channels, starting at 160, i.e. 160:191.
%    The numbering convention is that the first MEG channel is 0 (not 1).
if ~exist('triggerLinesPotential','var')
    triggerLinesPotential = [];
end
if isempty(triggerLinesPotential)
    triggerLinesPotential = 160:191;
end

% Create triggerHighThresh if not supplied or if equal to [].
if ~exist('triggerHighThresh','var')
    triggerHighThresh = [];
end
if isempty(triggerHighThresh)
    triggerHighThresh = 5; % mV
end

% Create triggerLowThresh if not supplied or if equal to [].
if ~exist('triggerLowThresh','var')
    triggerLowThresh = [];
end
if isempty(triggerLowThresh)
    triggerLowThresh = 0.2; % mV
end

triggerHighThresh_mV = triggerHighThresh*1000; % mV
triggerLowThresh_mV = triggerLowThresh*1000; % mV
clear triggerHighThresh triggerLowThresh

info = sqdread(sqdfile,'info');

fs = info.SampleRate; % Hz
triggerData = sqdread(sqdfile,'Channels',triggerLinesPotential);  % mV

% First find channels whose maximum is at least the "high" threshold and
% whose minimum is at most the "low" threshold.
triggerLinesTest = ( (max(triggerData) > triggerHighThresh_mV) & (min(triggerData) < triggerLowThresh_mV) );
triggerDataProbable = triggerData(:,triggerLinesTest);
triggerLinesProbable = triggerLinesPotential(triggerLinesTest);
clear triggerData

% (Not sure this is necessary.)
% Next, make sure that the time spent at rest is at least 100 times more
% than the time spent at the set value.
triggersNearHigh = mean(triggerDataProbable > triggerHighThresh_mV);
triggersNearLow = mean(triggerDataProbable < triggerLowThresh_mV);

triggerLinesHealthyHighRest =  (triggersNearHigh > 25 * triggersNearLow);
triggerDataHighRestGood = triggerDataProbable(:,triggerLinesHealthyHighRest);
triggerLinesHighRestGood = triggerLinesProbable(triggerLinesHealthyHighRest);
triggerLinesHighRestGoodCount = length(triggerLinesHighRestGood);

triggerLinesHealthyLowRest =  (triggersNearLow > 25 * triggersNearHigh);
triggerDataLowRestGood = triggerDataProbable(:,triggerLinesHealthyLowRest);
triggerLinesLowRestGood = triggerLinesProbable(triggerLinesHealthyLowRest);
triggerLinesLowRestGoodCount = length(triggerLinesLowRestGood);

if triggerLinesHighRestGoodCount > triggerLinesLowRestGoodCount
    triggerDataGood = triggerDataHighRestGood;
    triggerLinesGood = triggerLinesHighRestGood;
    triggerLinesGoodCount = triggerLinesHighRestGoodCount;
else
    triggerDataGood = triggerDataLowRestGood;
    triggerLinesGood = triggerLinesLowRestGood;
    triggerLinesGoodCount = triggerLinesLowRestGoodCount;
end


clear triggerDataProbable

% Inititialize outputs that hold trigger info
triggerCounts = nan(triggerLinesGoodCount,1);
triggersList = cell(triggerLinesGoodCount,1);

for iTrigger = 1:triggerLinesGoodCount
    if triggerLinesHighRestGoodCount > triggerLinesLowRestGoodCount
        triggerSet = find(triggerDataGood(:,iTrigger) < triggerLowThresh_mV);
    else
        triggerSet = find(triggerDataGood(:,iTrigger) > triggerHighThresh_mV);
    end
    % the next line uses find and diff to discard samples of triggerSet
    % that immediately follow any other values that cross
    % the "set" threshold.
    triggerSamples = (triggerSet([1; find(diff(triggerSet) > 1)+1]));
    triggerCounts(iTrigger) = length(triggerSamples);
    triggersList{iTrigger} = triggerSamples;
end

end




