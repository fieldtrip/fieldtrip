clear all
cbmex('open');
cbmex('trialconfig',1,'absolute')
  
samplingFreq = 30000;
eventChannel = 1;
contChannel = 1;
EventLFPs = [];
EventLFPTotal = [];
contRange = -50:50;
pause(0.5);

while size(EventLFPTotal,1) < 500
    [spikeData, procTime, contData] = cbmex('trialdata',1);
    %Example channel to choose for perievent marks
    continuousData = contData{contChannel,3};
    spikeTimes = spikeData{eventChannel,2}- procTime * samplingFreq;
    
    if ~isempty(spikeTimes)
        eventSpikes = bsxfun(@plus, contRange, double(spikeTimes(1:end-3)));
        EventLFPs = continuousData(eventSpikes);
        EventLFPTotal = [EventLFPTotal; EventLFPs];
        EventLFPTotal = reshape(EventLFPTotal, size(eventSpikes));        
        
    end
%    for i = 1:length(SpikeTimes)
        
%     	  try
%             EventLFPs(size(EventLFPs,1)+1,:) = ContinuousData(SpikeTimes(i):SpikeTimes(i)+100);
%         catch
%         end
%    end
    pause(0.5);
end
cbmex('close');

plot(mean(EventLFPTotal,1));
