function freq=fastpower(cfgf,data)
% This function is written for fast freq analysis based on
% freqanalysis_mtmconvol. It has the following fixed settings:
% cfgf.output       = 'pow';
% cfgf.keeptrials   = 'yes';
% cfgf.method       = 'mtmconvol';
% 
% The settings that need to be specified:
% cfgf.foi          = frequency of interest (recom: 10 )
% cfgf.toi          = time of interest      (recom: 0.5)
% cfgf.t_ftimwin    = timewindow            (recom: ones(1,length(cfgf.foi)) 
% cfgf.taper        = taper                 (recom: 'hanning')
%
% Copyright 2009, Ali Bahramisharif 

for i=1:length(data.time);
    data.offset(i) = round(data.time{i}(1)*data.fsample);
end
cfgf.trials = 'all';


% it is faster if you do not assign channel
if isfield(cfgf,'channel')
    cfgf.channel = channelselection(cfgf.channel, data.label);
else
    cfgf.channel=data.label;
end

sgnindx     = match_str(data.label, cfgf.channel);
numsgn      = size(sgnindx,1);
% if rectan is 1 it means that trials are of equal lengths
numper       = numel(data.trial);
numdatbnsarr = zeros(numper, 1);
for perlop = 1:numper
    numdatbnsarr(perlop) = size(data.trial{perlop},2);
end
min_smp = min(data.offset);
max_smp = max(numdatbnsarr(:)+data.offset(:));
cfgf.pad = (max_smp-min_smp) ./ data.fsample;
clear min_smp max_smp
numsmp = round(cfgf.pad .* data.fsample);
minoffset = min(data.offset);
timboi = round(cfgf.toi .* data.fsample - minoffset);
toi    = round(cfgf.toi .* data.fsample) ./ data.fsample;
numtoi = length(cfgf.toi);
numfoi = length(cfgf.foi);
numtap = zeros(numfoi,1);

% compute the tapers and their fft
knlspctrmstr = cell(numfoi,1);
for foilop = 1:numfoi
    acttapnumsmp = round(cfgf.t_ftimwin(foilop) .* data.fsample);
    % create a single taper according to the window specification as a replacement for the DPSS (Slepian) sequence
    tap = window(cfgf.taper, acttapnumsmp);
    tap = tap./norm(tap);
    % freqanalysis_mtmconvol always throws away the last taper of the Slepian sequence, so add a dummy taper
    tap(:,2) = nan;
    numtap(foilop) = size(tap,2)-1;
    ins = ceil(numsmp./2) - floor(acttapnumsmp./2);
    prezer = zeros(ins,1);
    pstzer = zeros(numsmp - ((ins-1) + acttapnumsmp)-1,1);
    ind    = (0:acttapnumsmp-1)' .* ((2.*pi./data.fsample) .* cfgf.foi(foilop));
    knlspctrmstr{foilop} = complex(zeros(numtap(foilop),numsmp));
    for taplop = 1:numtap(foilop)
        coswav  = vertcat(prezer,tap(:,taplop).*cos(ind),pstzer);
        sinwav  = vertcat(prezer,tap(:,taplop).*sin(ind),pstzer);
        wavelet = complex(coswav, sinwav);
        % store the fft of the complex wavelet
        knlspctrmstr{foilop}(taplop,:) = fft(wavelet,[],1)';
    end
end
powspctrm     = zeros(numper,numsgn,numfoi,numtoi);
numdatbns = numdatbnsarr(perlop,1);
prepad = zeros(numsgn,data.offset(perlop) - minoffset);
pstpad = zeros(numsgn,minoffset + numsmp - (data.offset(perlop) + numdatbns));
tmp = data.trial{perlop}(sgnindx,:);
tmp = [prepad tmp pstpad];
% avoid the use of a 3rd input argument to facilitate compatibility with star-P
% use explicit transpose, to avoid complex conjugate transpose
datspctra = transpose(fft(transpose(tmp)));

for foilop = 1:numfoi
    actfoinumsmp    = cfgf.t_ftimwin(foilop) .* data.fsample;
    acttimboiind    = find(timboi >= (-minoffset + data.offset(perlop) + (actfoinumsmp ./ 2)) & timboi <  (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
    nonacttimboiind = find(timboi <  (-minoffset + data.offset(perlop) + (actfoinumsmp ./ 2)) | timboi >= (-minoffset + data.offset(perlop) + numdatbns - (actfoinumsmp ./2)));
    acttimboi       = timboi(acttimboiind);
    numacttimboi    = length(acttimboi);
    for taplop = 1:numtap(foilop)
        autspctrmacttap = complex(zeros(numsgn,numacttimboi), zeros(numsgn,numacttimboi));
        if numacttimboi > 0
            for sgnlop = 1:numsgn
                dum = fftshift(ifft(datspctra(sgnlop,:) .* knlspctrmstr{foilop}(taplop,:),[],2));
                autspctrmacttap(sgnlop,:) = dum(acttimboi);
            end
        end
        
        powdum = 2.* abs(autspctrmacttap) .^ 2 ./ actfoinumsmp;
        powspctrm(perlop,:,foilop,acttimboiind) = powspctrm(perlop,:,foilop,acttimboiind) + reshape(powdum ./ numtap(foilop),[1,numsgn,1,numacttimboi]);
        powspctrm(perlop,:,foilop,nonacttimboiind) = nan;
    end
end

% collect the results
freq.label      = data.label(sgnindx);
freq.dimord     = 'rpt_chan_freq_time';
freq.freq       = cfgf.foi;
freq.time       = toi;
freq.powspctrm  = powspctrm;
freq.grad = data.grad;
freq.cfg = cfgf;
return