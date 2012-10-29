function [ orig ] = plx_orig_header( fname )

% PLX_ORIG_HEADER Extracts the header informations of plx files using the
% Plexon Offline SDK, which is available from
% http://www.plexon.com/assets/downloads/sdk/ReadingPLXandDDTfilesinMatlab-mexw.zip
%
% Use as
%   [orig] = plx_orig_header(filename)
%
% Copyright (C) 2012 by Thomas Hartmann
%
% This code can be redistributed under the terms of the GPL version 3 or
% newer.

% get counts...
[tscounts, wfcounts, evcounts, contcounts] = plx_info(fname, 0);

orig.TSCounts = tscounts;
orig.WFCounts = wfcounts;

% get event channels...
[dummy, evchans] = plx_event_chanmap(fname);

orig.EVCounts = zeros(512, 1);

for i=1:length(evchans)
  try
    orig.EVCounts(evchans(i)+1) = evcounts(i);
  catch
  end
end %for

% get more infos...
[OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = plx_information(fname);
orig.MagicNumber = 1480936528;
orig.Version = Version;
orig.ADFrequency = Freq;
orig.Comment = Comment;
orig.NumEventChannels = length(evchans);
orig.NumSlowChannels = length(contcounts);
[orig.NumDSPChannels, dummy] = plx_chanmap(fname);
orig.NumPointsWave = NPW;
orig.NumPointsPreThr = PreTresh;
[orig.Year orig.Month orig.Day orig.Hour orig.Minute orig.Second] = datevec(DateTime);
orig.LastTimestamp = Duration * Freq;
orig.Trodalness = Trodalness;
orig.DataTrodalness = Trodalness;
orig.BitsPerSpikeSample = SpikeADResBits;
orig.BitsPerSlowSample = SlowADResBits;
orig.SpikeMaxMagnitudeMV = SpikePeakV;
orig.SlowMaxMagnitudeMV = SlowPeakV;

% gather further info for additional headers...
[dummy, chan_names] = plx_chan_names(fname);
[dummy, dspchans] = plx_chanmap(fname);
[dummy, gains] = plx_chan_gains(fname);
[dummy, filters] = plx_chan_filters(fname);
[dummy, thresholds] = plx_chan_thresholds(fname);
[dummy, evnames] = plx_event_names(fname);
[dummy, ad_names] = plx_adchan_names(fname);
[dummy, adchans] = plx_ad_chanmap(fname);
[dummy, ad_freqs] = plx_adchan_freqs(fname);
[dummy, ad_gains] = plx_adchan_gains(fname);

% do channelheaders...
for i=1:orig.NumDSPChannels
  head = [];
  head.Name = ['mua_' chan_names(i, :)];
  head.SIGName = ['mua_' chan_names(i, :)];
  head.Channel = dspchans(i);
  head.SIG = i;
  head.Gain = gains(i);
  head.Filter = filters(i);
  head.Threshold = thresholds(i);
  
  orig.ChannelHeader(i) = head;
end %for

% do eventheaders...
for i=1:orig.NumEventChannels
  head = [];
  head.Name = evnames(i, :);
  head.Channel = evchans(i);
  
  orig.EventHeader(i) = head;
end %for

% do slowchannelheaders...
for i=1:orig.NumSlowChannels
  head = [];
  head.Name = ['lfp_' ad_names(i, :)];
  head.Channel = adchans(i);
  head.ADFreq = ad_freqs(i);
  head.Gain = ad_gains(i);
  
  orig.SlowChannelHeader(i) = head;
end %for

end
