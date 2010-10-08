% inital version of Matlab acquisition and saving tool for BIOSEMI amp
% TODO: add proper configuration, turn into function
% (C) 2010 S. Klanke 

try
	dummy = biosemix([0 0]);
catch me
	if strfind(me.message,'BIOSEMI device')
		clear biosemix
	else
		rethrow(me);
	end
end

S = biosemix;

numEEG = 32;
numAIB = 0;
cfg.ftOutput = 'buffer://localhost:1972';
cfg.gdfOutput = 'D:\biosemi.gdf';
cfg.decimate = 8;
cfg.order = 4;

numBlocks = 0;
numSamples = 0;



HDR.SampleRate = S.fSample;
HDR.NS = 1+numEEG+numAIB;
HDR.SPR = 1; 
HDR.NRec = 1; % always continuous
HDR.FileName = cfg.gdfOutput;
HDR.TYPE = 'GDF';
HDR.T0 = clock;

HDR.label   = cell(HDR.NS,1);
HDR.label{1} = 'Trigger';
for k=1:numEEG
	HDR.label{1+k} = sprintf('EEG%03i',k);
end
for k=1:numAIB
	HDR.label{1+numEEG+k} = sprintf('AIB%03i',k);
end

[datatyp,limits,datatypes,HDR.Bits,HDR.GDFTYP]=gdfdatatype('int32');

HDR.PhysDimCode = 512*ones(HDR.NS,1); % physicalunits('-')
HDR.DigMin = limits(1)*ones(HDR.NS,1);
HDR.DigMax = limits(1)*ones(HDR.NS,1);
HDR.PhysMin = HDR.DigMin;
HDR.PhysMax = HDR.DigMax;
HDR.FLAG.UCAL = 1;

HDR = sopen(HDR,'w');

[B,A] = butter(cfg.order, 0.8/cfg.decimate, 'low');
DM = online_downsample_init(cfg.decimate);
FM = online_filter_init(B, A, single(zeros(numEEG+numAIB,1)));


hdr.Fs = S.fSample;
hdr.nChans  = numEEG+numAIB;
hdr.label = HDR.Label(2:end);


ft_write_data(cfg.ftOutput, single(zeros(hdr.nChans,0)), 'header', hdr, 'append', false);
dummy = biosemix;	

while 1
	dat = biosemix([numEEG numAIB]);
	if isempty(dat)
		continue
	end
	N = size(dat,2);
	
	sig = single(dat(2:end,:)) * 0.262 / 2^31;
	
	[FM, fsig] = online_filter_apply(FM, sig);
	[DM, dsig] = online_downsample_apply(DM, fsig);
	
	ft_write_data(cfg.ftOutput, dsig, 'header', hdr, 'append', true);
	HDR = swrite(HDR, dat');
	flush;
	numBlocks = numBlocks + 1;
	numSamples = numSamples + N;
	fprintf(1,'%i blocks, %i samples\n', numBlocks, numSamples);
end