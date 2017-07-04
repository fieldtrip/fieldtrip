function [lcoh] = ft_connectivity_laggedcoherence(cfg,freqout)

% FT_CONNECTIVITY_LAGGEDCOHERENCE performs time-resolved coherence analysis 
% of oscillatory activity only, both within and between recording sites.
%
% Use as
%  lcoh=ft_connectivityanalysis(cfg,freqout) with cfg.method='laggedcoherence';
%  or as lcoh=ft_connectivity_laggedcoherence(cfg,freqout)
%
% The input data should be organised in a structure as obtained from
%  the FT_FREQANALYSIS function (freqout), such that freqout contains the
%  fields 'fourierspctrm' and 'time'.
%  The timepoints must be chosen such that the desired cfg.lag/cfg.foi
%  (lag in s) is an integer multiple of the time resolution in freqout.
%
%  At the moment, this function must be called separately for each frequency
%  of interest. To analyse multiple frequencies, we advise the use of a
%  forloop such as this:
%  cfg_F  = [];
%  cfg_LC = [];
%  cfg_F.method      = 'wavelet';
%  cfg_F.output      = 'fourier';
%  cfg_F.width       = 3;
%  cfg_F.keeptrials  = 'yes';
%  cfg_LC.lag        = cfg_F.width;
%  cfg_LC.method     = 'laggedcoherence';
%  foi               = 1:1:100;
%  fs                = data.fsample;
%  for counter = 1:length(foi);
%      cfg_F.foi     = foi(counter);
%      cfg_LC.foi    = foi(counter);
%      width         = cfg_F.width/cfg_F.foi;
%      cfg_F.toi     = data.time{1}(1) + ceil(fs*width/2)/fs : ...   %from:
%        cfg_LC.lag/cfg_F.foi : ...                     %in steps of size:
%        data.time{1}(end) - ceil(fs*width/2)/fs;                     %to:
%      freqout       = ft_freqanalysis(cfg_F,data);
%      lcoh(counter) = ft_connectivityanalysis(cfg_LC,freqout);
%  end
%
% The configuration for ft_connectivity_laggedcoherence should contain:
%   cfg.foi          =    frequency of interest (default=freqout.freq(1))
%   cfg.lag          =    the number of periods between the onset of the
%      time window used for phase estimate 1 and the onset of the time
%      window for phase estimate 2 (the default cfg.lag is set to match the
%      time resolution in freqout). We recommend users to choose cfg.lag such
%      that it is larger or equal to the width of the wavelet used for each
%      Fourier transform in ft_freqanalysis
%   cfg.output       =    'lcoh', or 'csd' (default='lcoh'). When the output
%      is set to 'csd', one can specify the channel combinations between
%      which to compute the lagged cross-spectra in cfg.channelcmb.
%
%   To calculate lagged coherence values from the cross-spectra, do:
%      abs(lcoh.laggedcrsspctrm)./sqrt(lcoh.powspctrm1.*lcoh.powspctrm2)
%      where lcoh.powspctrm1 denotes power in the channels of the first
%      column of cfg.channelcmb, and lcoh.powspctrm2 does the same for the
%      channels in the second column of cfg.channelcmb. Note that the power
%      is calculated for the same time windows that are used for calculating
%      the lagged cross-spectra.
%
% Optional settings:
%   cfg.timeresolved =    'yes' or 'no' (default='no'). If set to yes, lagged 
%      coherence is calculated separately for each pair of timepoints that
%      is separated by cfg.lag
%   cfg.nlags        =    The lags in lcoh are the set of (1:1:nlags)*lag
%      (default: cfg.nlag=1). Note that if cfg.timeresolved=='yes', then
%      cfg.nlags must be set to 1.
%   cfg.channel      =    Nx1 cell-array with selection of channels
%      (default = 'all'), see FT_CHANNELSELECTION for details
%   cfg.channelcmb   =   Mx2 cell-array with selection of channel pairs,
%      (default={'all','all'}), see FT_CHANNELCOMBINATION for details.
%   cfg.autocmb      =    'yes' or 'no' (default='no'). Adds all auto-
%      combinations of cfg.channel to cfg.channelcmb
%   cfg.trialsets    =    cell array with per cell the set of trials over
%      which lcoh is calculated. Each cell must contain 'all' or a 1xN
%      vector of trial indices. Default={'all'}. Note that this differs
%      from the required format of cfg.trials in e.g. ft_connectivityanalysis.
%
% When this measure is used for your publication, please cite:
% Fransen, Anne M. M, Van Ede, Freek, Maris, Eric (2015) Identifying
%  oscillations on the basis of rhythmicity. NeuroImage 118: 256-267.

% Copyright (C) 2015, Anne M.M. Fransen, DCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
% FieldTrip is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% FieldTrip is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',     {'timbinlag',         'lag'});
cfg = ft_checkconfig(cfg, 'renamed',     {'timbinlagincycles', 'lag'});
cfg = ft_checkconfig(cfg, 'renamed',     {'ntimbinlags',       'nlags'});

% ensure that the required options are present
cfg.channel             =   ft_getopt(cfg, 'channel',           'all');
cfg.trials              =   ft_getopt(cfg, 'trialsets',         {'all'});
cfg.output              =   ft_getopt(cfg, 'output',            'lcoh');
cfg.nlags               =   ft_getopt(cfg, 'nlags',             1);
cfg.foi                 =   ft_getopt(cfg, 'foi',               freqout.freq(1));
cfg.timeresolved        =   ft_getopt(cfg, 'timeresolved',      'no');
cfg.autocmb             =   ft_getopt(cfg, 'autocmb',           'no');
cfg.precision           =   ft_getopt(cfg, 'precision',          0.001); %timepoints that are separated by cfg.lag are identified with this precision (in s)
cfg.feedback            =   ft_getopt(cfg, 'feedback',          'no');
cfg.inputlock           =   ft_getopt(cfg, 'inputlock',         []);  % this can be used as mutex when doing distributed computation
cfg.outputlock          =   ft_getopt(cfg, 'outputlock',        []);  % this can be used as mutex when doing distributed computation
if strcmp(cfg.autocmb,'yes');                                   cfg.autocmb=true;
elseif strcmp(cfg.autocmb,'no');                                cfg.autocmb=false;
end
if ~iscell(cfg.trialsets); cfg.trialsets={cfg.trialsets}; end

% adjust cfg.foi to nearest possibility in freqout.freq, and only allow 1 foi
[dummy indexf]          =   min(abs(freqout.freq - cfg.foi(1)));
cfg.foi                 =   freqout.freq(indexf);

% check if the input data is valid for this function
freqout = ft_checkdata(freqout, 'datatype', {'freq'}, 'feedback', cfg.feedback,'dimord',{'rpttap_chan_freq_time','rpt_chan_freq_time'});
if ~isfield(freqout,'fourierspctrm')
    ft_error('ft_connectivity_laggedcoherence requires frequency data with a fourierspctrm.');
elseif size(freqout.fourierspctrm,4)==1
    ft_error('ft_connectivity_laggedcoherence requires frequency data with a time axis')
end

% ensure that cfg.lag is present
cfg.lag  =  ft_getopt(cfg, 'lag',  (freqout.time(2:2:end)-freqout.time(1:2:end-1))*cfg.foi);

% select channels of interest
tmpcfg = [];
tmpcfg.channel          =  cfg.channel;
freqout                 =  ft_selectdata(tmpcfg, freqout);
cfg.channel             =  cell(length(freqout.label),1);
cfg.channel(:)          =  freqout.label;

% select trials of interest for each trialset
for counter=1:length(cfg.trialsets)
    if strcmp(cfg.trialsets{counter},'all')
        cfg.trialsets{counter}   =    1:size(freqout.fourierspctrm,1);
    end
    if size(cfg.trialsets{counter},1)==1
        cfg.trialsets{counter}=cfg.trialsets{counter}';
    end
end

% some proper error handling
if numel(freqout.label)==0
    ft_error('no channels were selected');
end
if ~(strcmp(cfg.output,'lcoh') || strcmp(cfg.output,'csd'))
    ft_error('cfg.output must be either ''lcoh'' or ''csd''');
end
index = cellfun('isclass',cfg.trialsets,'char');
if sum(index)>0
    index2               =  cellfun(@(x) strcmp(x,'all'),cfg.trialsets(index),'UniformOutput',true);
    if any(index2==0)
        ft_error('each cell of cfg.trialsets must contain either a 1xN vector, or ''all''');
    end
    cfg.trialsets(index) =  {1:size(freqout.fourierspctrm,1)};
end
clear tmpcfg index index2

% prepare cfg.channelcmb
csdflag = strcmp(cfg.output,'csd');
if ~csdflag
    %only autocombinations
    cfg.channelcmb      = [cfg.channel,cfg.channel];
else
    if ~isfield(cfg, 'channelcmb')
        %set the default for the channelcombination
        cfg.channelcmb  = {'all' 'all'};
    end
    cfg.channelcmb      = ft_channelcombination(cfg.channelcmb,freqout.label,cfg.autocmb);
end
% determine the corresponding indices of all channel combinations
[dummy,chancmbind(:,1)] = match_str(cfg.channelcmb(:,1), cfg.channel);
[dummy,chancmbind(:,2)] = match_str(cfg.channelcmb(:,2), cfg.channel);
clear dummy
if csdflag
    assert(length(unique(chancmbind(:)))>1, 'CSD output requires multiple channels');
end

%initiate some variables
nchancmb    =   size(chancmbind,1);
nrep        =   size(freqout.fourierspctrm,1);
ntrialsets  =   length(cfg.trialsets);
cyclelength =   1/cfg.foi;
dtim        =   repmat(freqout.time,fliplr(size(freqout.time))) - repmat(freqout.time',size(freqout.time));
dtim        =   single(round(dtim./cfg.precision).*cfg.precision);

% switch between time-resolved and non-time-resolved calculation of lagged coherence
switch cfg.timeresolved
    case 'no'
        %initiate some variables
        laggedcps       = complex(zeros(nchancmb,cfg.nlags,nrep));
        power           = complex(zeros(nchancmb,cfg.nlags,2,nrep));
        hasdata         = false(nrep,cfg.nlags);
        nsmplslaggedcps = zeros(cfg.nlags,nrep);
        lagwidth        = zeros(1,cfg.nlags);
        
        for lagindx=1:cfg.nlags
            % select all pairs of timepoints with relative lag cfg.lag (identified with precision cfg.precision)
            lagwidth(lagindx) = cfg.lag*cyclelength*lagindx;
            lagwidth(lagindx) = dtim(find(abs(dtim-lagwidth(lagindx))==min(abs(dtim(:)-lagwidth(lagindx))),1));
            [t1 t2]           = find(dtim==lagwidth(lagindx));
            
            % calculate laggedcrossproducts and power per channelcmb and trial
            for trialindx=1:nrep
                % get the spectrum for this trial and frequency
                fcs1         = complex(zeros(size(freqout.fourierspctrm,2),length(t1)));
                fcs2         = fcs1;  % fcs1 and fcs2 have dimord chan_time
                fcs1(:,:)    = freqout.fourierspctrm(trialindx,:,indexf,t1);
                fcs2(:,:)    = freqout.fourierspctrm(trialindx,:,indexf,t2);
                colswithnans = any(isnan(fcs1)|isnan(fcs2),1);
                fcs1(:,colswithnans) = []; fcs2(:,colswithnans) = [];
                % sum laggedcrossproducts and power over all timepoints
                for tcounter=1:length(t1)
                    laggedcrossproduct = fcs1(chancmbind(:,1),tcounter).*conj(fcs2(chancmbind(:,2),tcounter));
                    laggedcps(:,lagindx,trialindx)     = laggedcps(:,lagindx,trialindx)+ laggedcrossproduct;
                    power(:,lagindx,1,trialindx)       = power(:,lagindx,1,trialindx)+ abs(fcs1(chancmbind(:,1),tcounter)).^2;
                    power(:,lagindx,2,trialindx)       = power(:,lagindx,2,trialindx)+ abs(fcs2(chancmbind(:,2),tcounter)).^2;
                    hasdata(trialindx,lagindx)         = true;
                    nsmplslaggedcps(lagindx,trialindx) = nsmplslaggedcps(lagindx,trialindx)+1;
                end
            end
        end
        
        % calculate lagged coherence
        if strcmp('lcoh',cfg.output)
            laggedcoh         = complex(nan(ntrialsets,nchancmb,cfg.nlags));
            for trialset=1:ntrialsets
                % sum laggedcps and laggedpower per trialset
                trials        = cfg.trialsets{trialset};
                sumlaggedcps  = nansum(laggedcps(:,:,trials),3);
                sumpower      = nansum(power(:,:,:,trials),4);
                % normalise the lagged csd
                denomcoh      = sqrt(sumpower(:,:,1).*sumpower(:,:,2));
                laggedcoh(trialset,:,:) = abs(sumlaggedcps)./denomcoh;
                for lagindex=1:cfg.nlags
                    laggedcoh(trialset,:,sum(hasdata(trials,lagindx),1)==0) = NaN;
                end
            end
            % make output structure
            lcoh = [];
            lcoh.laggedcoh          = laggedcoh;
            lcoh.label              = freqout.label;
            lcoh.freq               = cfg.foi;
            lcoh.lag                = lagwidth*cfg.foi;
            if length(cfg.trialsets)==1
                lcoh.laggedcoh      = shiftdim(lcoh.laggedcoh,1);
                if length(lcoh.lag)==1
                    lcoh.dimord     = 'chan';
                else
                    lcoh.dimord     = 'chan_rep';
                end
            else
                if all(cellfun(@length,cfg.trialsets)==1)
                    trials=cell2mat(cfg.trialsets);
                    if isfield(freqout,'trialinfo')
                        lcoh.trialinfo=freqout.trialinfo(trials,:);
                    end
                else
                    if isfield(freqout,'trialinfo')
                        for counter=1:ntrialsets
                            lcoh.trialinfos{counter}=freqout.trialinfo(cfg.trialsets{counter},:);
                        end
                    end
                    lcoh.trialsets   = cfg.trialsets;
                end
                if length(lcoh.lag)==1
                    lcoh.dimord      = 'rep_chan';
                else
                    lcoh.dimord      = 'rep_chan_rep';
                end
            end
            
        % calculate lagged crossspectra and corresponding power spectra
        elseif csdflag
            laggedcrsspctrm    = complex(zeros(ntrialsets,nchancmb,cfg.nlags));
            powspctrm1         = complex(zeros(ntrialsets,nchancmb,cfg.nlags));
            powspctrm2         = complex(zeros(ntrialsets,nchancmb,cfg.nlags));
            for trialset=1:ntrialsets
                % sum laggedcps and laggedpower per trialset
                trials         = cfg.trialsets{trialset};
                sumlaggedcps   = nansum(laggedcps(:,:,trials),3);
                sumlaggedpower = nansum(power(:,:,:,trials),4);
                for lagindx=1:cfg.nlags
                    sumnsmpls  = nansum(nsmplslaggedcps(lagindx,trials),2);
                    laggedcrsspctrm(trialset,:,lagindx) = sumlaggedcps(:,lagindx)/sumnsmpls;
                    powspctrm1(trialset,:,lagindx)      = sumlaggedpower(:,lagindx,1)/sumnsmpls;
                    powspctrm2(trialset,:,lagindx)      = sumlaggedpower(:,lagindx,2)/sumnsmpls;
                    if sum(hasdata(trials,lagindx),1)==0
                        laggedcrsspctrm(trialset,:,lagindx) = NaN;
                    end
                end
            end
            % make output structure
            lcoh = [];
            lcoh.laggedcrsspctrm  = laggedcrsspctrm;
            lcoh.powspctrm1       = powspctrm1;
            lcoh.powspctrm2       = powspctrm2;
            lcoh.label            = freqout.label;
            lcoh.freq             = cfg.foi;
            lcoh.lag              = lagwidth*cfg.foi;
            lcoh.labelcmb         = cfg.channelcmb;
            if length(cfg.trialsets)==1
                lcoh.laggedcrsspctrm = shiftdim(lcoh.laggedcrsspctrm,1);
                lcoh.powspctrm1   = shiftdim(lcoh.powspctrm1,1);
                lcoh.powspctrm2   = shiftdim(lcoh.powspctrm2,1);
                if length(lcoh.lag)==1
                    lcoh.dimord   = 'chan';
                else
                    lcoh.dimord   = 'chan_rep';
                end
            else
                if all(cellfun(@length,cfg.trialsets)==1)
                    trials=cell2mat(cfg.trialsets);
                    if isfield(freqout,'trialinfo')
                        lcoh.trialinfo=freqout.trialinfo(trials,:);
                    end
                else
                    if isfield(freqout,'trialinfo')
                        for counter=1:ntrialsets
                            lcoh.trialinfos{counter}=freqout.trialinfo(cfg.trialsets{counter},:);
                        end
                    end
                    lcoh.trialsets = cfg.trialsets;
                end
                if length(lcoh.lag)==1
                    lcoh.dimord    = 'rep_chan';
                else
                    lcoh.dimord    = 'rep_chan_rep';
                end
            end
        end
        
    case 'yes'
        % method-specfic checks
        if cfg.nlags>1
            ft_error('when calculating timeresolved lcoh, cfg.nlags must be set to 1');
        end
        
        % select pairs of timepoints with relative lag cfg.lag (identified with precision cfg.precision)
        lagwidth            =  cfg.lag*cyclelength;
        lagwidth            =  dtim(find(abs(dtim-lagwidth)==min(abs(dtim(:)-lagwidth)),1));
        [t1 t2]             =  find(dtim==lagwidth);
        
        % initiate some variables (dependent on nr of timepoints)
        ntoi                =   length(t1);
        laggedcrossproduct  =   complex(zeros(nchancmb,ntoi,nrep));
        power               =   complex(zeros(nchancmb,ntoi,2,nrep));
        hasdata             =   false(nrep,ntoi);
        
        % calculate laggedcrossproducts and power per channel, time of interest, and trial
        for trialindx=1:nrep
            % get the spectrum for this trial and frequency
            fcs1 = complex(zeros(size(freqout.fourierspctrm,2),length(t1))); 
            fcs2 = fcs1;  % fcs1 and fcs2 have dimord chan_time
            fcs1(:,:)       = freqout.fourierspctrm(trialindx,:,indexf,t1);
            fcs2(:,:)       = freqout.fourierspctrm(trialindx,:,indexf,t2);
            colswithnans    = any(isnan(fcs1)|isnan(fcs2),1);
            fcs1(:,colswithnans) = [];
            fcs2(:,colswithnans) = [];
            % sum laggedcps and laggedpower over all timepoints
            for tcounter=1:ntoi
                laggedcrossproduct(:,tcounter,trialindx)= fcs1(chancmbind(:,1),tcounter).*conj(fcs2(chancmbind(:,2),tcounter));
                power(:,tcounter,1,trialindx)           = abs(fcs1(chancmbind(:,1),tcounter)).^2;
                power(:,tcounter,2,trialindx)           = abs(fcs2(chancmbind(:,2),tcounter)).^2;
                hasdata(trialindx,tcounter)             = true;
            end
        end
        
        % calculate lagged coherence
        if strcmp('lcoh',cfg.output)
            ntrialsets  = length(cfg.trialsets);
            laggedcoh   = complex(nan(ntrialsets,nchancmb,ntoi));
            for trialset=1:ntrialsets
                % sum laggedcps and laggedpower per trialset
                trials = cfg.trialsets{trialset};
                trialsumlaggedcps = nansum(laggedcrossproduct(:,:,trials),3);
                trialsumpower     = nansum(power(:,:,:,trials),4);
                for toii=1:ntoi
                    sumlaggedcps  = trialsumlaggedcps(:,toii);
                    sumpower      = trialsumpower(:,toii,:);
                    % normalise the power per frequency
                    denomcoh      = sqrt(sumpower(:,:,1).*sumpower(:,:,2));
                    laggedcoh(trialset,:,toii) = abs(sumlaggedcps)./denomcoh;
                end
            end
            % make output structure
            lcoh = [];
            lcoh.laggedcoh     = laggedcoh;
            lcoh.label         = freqout.label;
            lcoh.freq          = cfg.foi;
            lcoh.time          = nanmean([freqout.time(t1);freqout.time(t2)],1);
            lcoh.lag           = (1:1:cfg.nlags)*lagwidth*cfg.foi;
            if length(cfg.trialsets)==1
                lcoh.laggedcoh = shiftdim(lcoh.laggedcoh,1);
                lcoh.dimord    = 'chan_time';
            else
                if all(cellfun(@length,cfg.trialsets)==1)
                    trials=cell2mat(cfg.trialsets);
                    if isfield(freqout,'trialinfo')
                        lcoh.trialinfo=freqout.trialinfo(trials,:);
                    end
                else
                    if isfield(freqout,'trialinfo')
                        for counter=1:ntrialsets
                            lcoh.trialinfos{counter} = freqout.trialinfo(cfg.trialsets{counter},:);
                        end
                    end
                    lcoh.trialsets = cfg.trialsets;
                end
                lcoh.dimord      = 'rep_chan_time';
            end
            
        % calculate lagged crossspectra and corresponding power spectra
        elseif csdflag
            laggedcrsspctrm      = complex(nan(ntrialsets,nchancmb,ntoi));
            powspctrm1           = complex(nan(ntrialsets,nchancmb,ntoi));
            powspctrm2           = complex(nan(ntrialsets,nchancmb,ntoi));
            for trialset=1:ntrialsets
                % sum laggedcps and laggedpower per trialset
                trials           = cfg.trialsets{trialset};
                sumlaggedcps     = nansum(laggedcrossproduct(:,:,trials),3);
                sumpower         = nansum(power(:,:,:,trials),4);
                sumnsmpls        = length(trials);
                for toii=1:ntoi
                    laggedcrsspctrm(trialset,:,toii)= sumlaggedcps(:,toii)/sumnsmpls;
                    powspctrm1(trialset,:,toii)     = sumpower(:,toii,1)/sumnsmpls;
                    powspctrm2(trialset,:,toii)     = sumpower(:,toii,2)/sumnsmpls;
                end
            end
            % make output structure
            lcoh = [];
            lcoh.laggedcrsspctrm = laggedcrsspctrm;
            lcoh.powspctrm1      = powspctrm1;
            lcoh.powspctrm2      = powspctrm2;
            lcoh.label           = freqout.label;
            lcoh.freq            = cfg.foi;
            lcoh.time            = nanmean([freqout.time(t1);freqout.time(t2)],1);
            lcoh.lag             = (1:1:cfg.nlags)*lagwidth*cfg.foi;
            lcoh.labelcmb        = cfg.channelcmb;
            if length(cfg.trialsets)==1
                lcoh.laggedcrsspctrm = shiftdim(lcoh.laggedcrsspctrm,1);
                lcoh.powspctrm1  = shiftdim(lcoh.powspctrm1,1);
                lcoh.powspctrm2  = shiftdim(lcoh.powspctrm2,1);
                lcoh.dimord      = 'chan_time';
            else
                if all(cellfun(@length,cfg.trialsets)==1)
                    trials=cell2mat(cfg.trialsets);
                    if isfield(freqout,'trialinfo')
                        lcoh.trialinfo=freqout.trialinfo(trials,:);
                    end
                else
                    if isfield(freqout,'trialinfo')
                        for counter=1:ntrialsets
                            lcoh.trialinfos{counter} = freqout.trialinfo(cfg.trialsets{counter},:);
                        end
                    end
                    lcoh.trialsets = cfg.trialsets;
                end
                lcoh.dimord    = 'rep_chan_time';
            end
        end
end
return

