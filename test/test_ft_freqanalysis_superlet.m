function test_ft_freqanalysis_superlet

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis

% Generates toy data containing oscillation packets at 20 Hz, 40 Hz and 60 Hz,
% each with 2 "neighbour packets" at a neighbouring frequency (+10 Hz) and a 
% neighbouring time (+ 12 cycles).
% Performs additive superlet analysis on toy data using frequency-adaptive 
% superlet orders from 1 to 30 with base wavelet width 3 cycles.
% Reproduces Figure 3d1 in https://doi.org/10.1101/583732 

%% generate toy data 
% The code used to generate the toy data is taken from
% http://muresanlab.tins.ro/sources/superlets/superlets.zip
% with permission by the author

%general settings
fs = 1000;
N = fs * 3.5;

%generate channel frequency and time neighbors 
vfTarget =                  [20, 40, 60];   %target frequencies
vfNeighbF=                  [30, 50, 70];   %neighboring frequencies
nFreqs =                    numel(vfTarget);
nWaveCycles =               11;
nTNeighbSpacingInBursts=    1 + 2/nWaveCycles;
bTNeighbRelativeSpacing =   true;
nPacketSpacinginBursts =    1/nWaveCycles;
nPackets =                  2;
nPreSpaceS =                0.25;

nLongestBurstlen =  round(nWaveCycles * fs / min([vfTarget vfNeighbF]));
nPacketSpacing =    round(nPacketSpacinginBursts * nLongestBurstlen);
nPacketLen =        round(nLongestBurstlen * (2 + nTNeighbSpacingInBursts - 1));

xTrg = zeros(1,N);  %a place to store the target
xNF = zeros(1,N);   %a place to store the neighbor in frequency
xNT = zeros(1,N);   %a place to store the neighbor in time


for i = 1 : numel(vfTarget)         %for al target frequenc
    nPacketOffset = nPreSpaceS * fs + (i - 1)*(nPacketLen + nPacketSpacing);    %initial offset - useful to place at different time moments packets  
    fTarg =  vfTarget(i);
    fNeigF = vfNeighbF(i);
    nTarg =  round(nWaveCycles * (fs / fTarg));   %length of target samples
    vTarg =  sin(2*pi*fTarg/fs  * (0 : nTarg-1)); %target burst

    vfNeighF = linspace(fNeigF, fTarg, nPackets);   %neighboring frequencies
    if bTNeighbRelativeSpacing %same spacing no matter the frequency
        vfNeighT = linspace(round(nTNeighbSpacingInBursts * nTarg),0, nPackets); %dictated by the lowest frequency
    end

    for p = 1 : nPackets - 1
        fNeigF =    vfNeighF(p);
        nNeigF =    nTarg;    %length of the freq neighbor burst   
        vNeighbF =  sin(2*pi*fNeigF/fs * (0 : nNeigF-1) - pi/1.5);

        nFirstTargSample =  nPacketOffset +  (p - 1) * nFreqs * (nPacketLen + nPacketSpacing) + round((nLongestBurstlen-nTarg ) / 2);
        nFirstNFSample =    nFirstTargSample;
        nFirstNTSample =    nFirstTargSample + round(vfNeighT(p));
        
        xTrg(nFirstTargSample : nFirstTargSample + nTarg -  1) = vTarg; 
        xNF (nFirstNFSample   : nFirstNFSample +   nNeigF - 1) = vNeighbF;
        xNT (nFirstNTSample   : nFirstNTSample  +  nTarg -  1) = vTarg;
    end
end

xSignal = xTrg + xNF + xNT;

%% put data into a struct that Fieldtrip can work with
data.trial{1} = xSignal;
data.time{1} = 0:1/fs:3.5-1/fs;
data.label{1} = 'chan';

%% do spectral analysis
cfg = [];
cfg.method = 'superlet';
cfg.output = 'pow';
cfg.pad = 'nextpow2';
cfg.foi = 10:0.5:80;
cfg.toi = 0:0.001:3.5-1/fs;
cfg.superlet.basewidth = 3;
cfg.superlet.gwidth = 3;
cfg.superlet.combine = 'additive';
cfg.superlet.order = round(linspace(1,30,numel(cfg.foi)));

pow = ft_freqanalysis(cfg, data);

%% plot
figure()
imagesc(pow.time, pow.freq, squeeze(pow.powspctrm));
set(gca, 'YDir', 'normal');
xlabel('time [s]');
ylabel('freq [Hz]');
