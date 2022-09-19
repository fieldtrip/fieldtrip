function [NEV] = mergeNEV()

%% 
% Saves a new NEV file that contains event data from one NEV and spike data
% from another
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use: mergeNEV()

%mergeNEV version = '1.0.0.0';

%Author: Nick Halper
%Contact: nhalper@blackrockmicro.com

%% 
% Choose file that contains event data


uiwait(msgbox('Choose the file containing event data (comments, digital inputs, etc) that you want to keep','Choose Event Data','modal'));
EventNEV = openNEV();

uiwait(msgbox('Choose the file containing sorted spike data that you want to keep','Choose Spike Data','modal'));
SpikeNEV = openNEV();




for i = 1:length(SpikeNEV.MetaTags.ChannelID)
    EventIndices = find(EventNEV.Data.Spikes.Electrode == SpikeNEV.MetaTags.ChannelID(i));
    SpikeDataIndices = find(SpikeNEV.Data.Spikes.Electrode == SpikeNEV.MetaTags.ChannelID(i));
    
    EventNEV.Data.Spikes.TimeStamp(EventIndices) = [];
    EventNEV.Data.Spikes.Electrode(EventIndices) = [];
    EventNEV.Data.Spikes.Unit(EventIndices) = [];
    EventNEV.Data.Spikes.Waveform(:,EventIndices) = [];
    
    
    
    disp('Length of Spike Indices:');
    disp(length(SpikeDataIndices));
    
    disp('Length of Event Indices:');
    disp(length(EventIndices));
    
    EventNEV.Data.Spikes.TimeStamp = [EventNEV.Data.Spikes.TimeStamp SpikeNEV.Data.Spikes.TimeStamp(SpikeDataIndices)];
    EventNEV.Data.Spikes.Electrode = [EventNEV.Data.Spikes.Electrode SpikeNEV.Data.Spikes.Electrode(SpikeDataIndices)];
    EventNEV.Data.Spikes.Unit = [EventNEV.Data.Spikes.Unit SpikeNEV.Data.Spikes.Unit(SpikeDataIndices)];
    EventNEV.Data.Spikes.Waveform = [EventNEV.Data.Spikes.Waveform SpikeNEV.Data.Spikes.Waveform(:,SpikeDataIndices)];
    
    [EventNEV.Data.Spikes.TimeStamp, ISort] = sort(EventNEV.Data.Spikes.TimeStamp);
    EventNEV.Data.Spikes.Electrode = EventNEV.Data.Spikes.Electrode(ISort);
    EventNEV.Data.Spikes.Unit = EventNEV.Data.Spikes.Unit(ISort);
    EventNEV.Data.Spikes.Waveform = EventNEV.Data.Spikes.Waveform(:,ISort);
    
end


SpikeNEV.MetaTags.ChannelID = EventNEV.MetaTags.ChannelID;
SpikeNEV.MetaTags.HeaderOffset = EventNEV.MetaTags.HeaderOffset;
EventNEV.MetaTags = SpikeNEV.MetaTags;


%EventNEV.Data.Spikes = SpikeNEV.Data.Spikes;

NEV = EventNEV;


uiwait(msgbox('Choose the save location for the new file','Choose Save Location','modal'));
[FileName,PathName] = uiputfile('.nev');


SavePath = fullfile(PathName,FileName);

saveNEV(NEV,fullfile(PathName,FileName));

end

