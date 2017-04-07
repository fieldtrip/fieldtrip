function test_bug905

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_freqanalysis ft_specest_mtmfft

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug905.mat'));

timeunitinsec=10^(-4); % in seconds
fs=1/(4.8*timeunitinsec);
timestep=1/fs;

windowlength=2;
windowlengthinsamples=ceil(windowlength/timestep);
windowlengthinsec=timestep*windowlengthinsamples;

% determine the fourier frequencies
ralfreq=1/windowlengthinsec;
minfreq=ralfreq;
maxfreq=150;
freqs=minfreq:ralfreq:maxfreq;
nfreq=length(freqs);

smoothfraction=0.50;
smoothfreqs=freqs*smoothfraction;
ntaps=floor(smoothfreqs*windowlengthinsec - 1);
% identify the first frequency bin
selvec=(ntaps<=1);
freqsthisbin=freqs(selvec);
freqbins{1}=[freqsthisbin(1),freqsthisbin(end)];
% identify the remaining frequency bins
% we make frequency bins of increasing width
nfreqsinfreqbin=length(freqsthisbin);
lastfreqindxthisbin=nfreqsinfreqbin;
finished=false;
freqbinindx=1;
while ~finished
    freqbinindx=freqbinindx+1;
    nfreqsinfreqbin=nfreqsinfreqbin+2;
    firstfreqindxthisbin=lastfreqindxthisbin+1;
    lastfreqindxthisbin=lastfreqindxthisbin+nfreqsinfreqbin;
    if lastfreqindxthisbin>=nfreq
        lastfreqindxthisbin=nfreq;
        finished=true;
    end;
    freqsthisbin=freqs(firstfreqindxthisbin:lastfreqindxthisbin);
    freqbins{freqbinindx}=[freqsthisbin(1),freqsthisbin(end)];
end;
nfreqbins=length(freqbins);
tapsmofreqs=zeros(1,nfreqbins);
freqbinindcs=zeros(1,nfreq);
timeresol=1/fs;
mintapsmofreq=1/(windowlengthinsec - timeresol);
for freqbinindx=1:nfreqbins
    tmpval=mean(freqbins{freqbinindx})*smoothfraction/2;
    tapsmofreqs(freqbinindx)=max([tmpval,mintapsmofreq]);
    freqindx1=find(freqs==freqbins{freqbinindx}(1));
    freqindx2=find(freqs==freqbins{freqbinindx}(2));
    freqbinindcs(freqindx1:freqindx2)=freqbinindx;
end;

% perform a frequency analysis
allfreqs=[];
freqout_mtmfft=cell(1,nfreqbins);
for freqbinindx=1:nfreqbins
  freqcfg=[];
  freqcfg.method = 'mtmfft';
  freqcfg.output = 'fourier';
  freqcfg.keeptrials = 'yes';
  freqcfg.foilim = freqbins{freqbinindx};
  freqcfg.tapsmofrq = tapsmofreqs(freqbinindx);
  freqcfg.taper = 'dpss';
  freqcfg.pad='maxperlen';
  freqout_mtmfft{freqbinindx}=ft_freqanalysis(freqcfg,datapart);
  allfreqs=[allfreqs,freqout_mtmfft{freqbinindx}.freq];
end;

ntrials=length(datapart.trial);
trialindcs=1:ntrials;

freqdescrout_mtmfft=cell(1,nfreqbins);
for freqbinindx=1:nfreqbins
  freqdescrcfg=[];
  freqdescrcfg.trials=trialindcs;
  freqdescrout_mtmfft{freqbinindx}=ft_freqdescriptives(freqdescrcfg,freqout_mtmfft{freqbinindx});
end;


