function test_old_ft_sourceanalysis_lcmv

% MEM 1gb
% WALLTIME 00:10:00

% TEST test_old_ft_sourceanalysis_lcmv

%% Test script, by Johanna

ft_defaults;
swd='/data/johzum/delange2008data/';
% this resides on both mentat248 and mentat269
% the raw data is from Biomag 2010 competition
cd(swd);
cfg.dataset = '/data/johzum/delange2008data/dataset02.ds';

% specify trial lengths/types
cfg.trialfun='trialfun_prestim_postresp';
cfg.trialdef.preeventtype  = 'backpanel trigger';
cfg.trialdef.preeventvalue = [41 42 43 44 45 51 52 53 54 55];
cfg.trialdef.prestim    = 1; % seconds prior to stimulus onset
cfg.trialdef.posteventtype  = 'backpanel trigger';
cfg.trialdef.posteventvalue = [1 2 11 12];
cfg.trialdef.postresp=1; % seconds after subject response
cfg.trialdef.stimrespmaxgap=4.0;
[cfg] = ft_definetrial(cfg)
% data.trialinfo 2nd column is now resptype


% load data into trials
cfg.padding=0; %total length greater than longest trial, in seconds
cfg.continuous='yes';
[data] = ft_preprocessing(cfg);
save(['/data/johzum/delange2008data/cfg2.mat'],'cfg');
save(['/data/johzum/delange2008data/datastim.mat'],'-v7.3','-struct','data');
% rt is now 1st column of data.trialinfo
% for ii=1:length(data.time),endpt(ii)=data.time{ii}(end);end
% rt=endpt-.8;
% [sortrt,rtind]=sort(rt);

% create lead field
hdr=ft_read_header(cfg.dataset);
cfg.grad=hdr.grad;
cfg.headshape='dataset02.shape';
ft_headshape=ft_prepare_localspheres(cfg);
cfg.headshape=[];
cfg.vol=ft_headshape;
save(['/data/johzum/delange2008data/ft_headshape.mat'],'ft_headshape');
cfg.grid.resolution=.8;
grid=ft_prepare_leadfield(cfg,data);
save(['/data/johzum/delange2008data/grid.mat'],'grid');
save(['/data/johzum/delange2008data/cfg4.mat'],'cfg');


%% nutmeg-style time-freq lcmv beamformer
clear all global
load grid.mat
load ft_headshape.mat

try % if already run through before
    load('leadfield1.mat')
catch % if going through first time
    cfg=[];
    cfg.vol=ft_headshape;
    cfg.dataset = '/data/johzum/delange2008data/dataset02.ds';
    hdr=ft_read_header(cfg.dataset);
    [megallchan] = ft_channelselection({'MEG' 'MEGREF'},hdr.label);
    cfg.channel=megallchan;
    [vol,grad]=ft_prepare_vol_sens(cfg.vol,hdr.grad,'channel',hdr.label);
    % leadfield = ft_compute_leadfield(grid.pos(grid.inside,:), grad, vol, 'reducerank', 'no', 'normalize', 'no');
    leadfield1 = ft_compute_leadfield(grid.pos, grad, vol, 'reducerank', 'no', 'normalize', 'no');
    % save('leadfield.mat','leadfield');
    save('leadfield1.mat','leadfield1');
end

ctf_ss=ft_read_vol('dataset02.hdm');
cfg=[];
cfg.method='fiducial';
cfg.fiducial=ctf_ss.mri;
% this analyze format below was created from the .mri to which the fiducials
% match. first save-as the .mri to V2.  then use older CTF software save to
% analyze. see:
% http://nutmeg.berkeley.edu/index.php?title=Getting_Started#Use_CTF_.mri_to_which_.2A.hdm_file_already_contains_fiducial_information
mriav2=ft_read_mri('/data/johzum/delange2008data/dataset02_V2.img');
mriav2r=ft_volumerealign(cfg,mriav2);

cfg.grid=grid;
lfrs=reshape(leadfield1,[size(leadfield1,1) 3 size(cfg.grid.pos,1)]);
for ii=1:size(lfrs,3)
    grid.leadfield{ii}=lfrs(:,:,ii);
end

freq=[8 12; 12 25; 25 48; 52 90; 90 130];
if ~exist('data','var')
    data=load('datastim.mat');
end
for jj=1:5
    cfg=[];
    cfg.padding=0; %not needed after initial load.
    cfg.bpfilter='yes';
    cfg.demean='yes';
    %%% I would prefer to use Nutmeg's call to firls filter via nut_filter2.m, 
    %%% even though it uses signal toolbox.
    cfg.baselinewindow=[-.8 -0.1]; % exclude first 200ms from demean in case of filter artefacts
    cfg.bpfiltord=3;
    cfg.bpfiltertype='but'; %add nutmeg firls?
    cfg.bpfreq=freq(jj,:);
    cfg.keeptrials='yes';
    cfg.channel={'MEG' 'MEGREF'};
    % trial 629 wonky for channel 102 during prestim.  throw out trial 629.
    % channel 96 is wonk also, trial 741
    cfg.trials=[1:628 630:740 742:799];
    [freqdata] = ft_preprocessing(cfg, data);
    save(['freqdata' num2str(jj) '.mat'],'-v7.3','freqdata');
    clear freqdata
end
clear data

%% stim-locked analysis
for jj=1:5
    flags.avetime=1;
    load(['freqdata' num2str(jj) '.mat'])
    winstim=[freqdata.time{1}(1)+[.2 .4:.1:2.1]; freqdata.time{1}(1)+[.6 .8:.1:2.5]]';
    timefreqstim{size(winstim,1),size(winstim,1)}=[];
    for kk=[2:size(winstim,1)]
        cfg=[];
        cfg.covariance='yes';
        cfg.permutation='no';
        cfg.vartrllength=1;
        cfg.keeptrials='no';
        cfg.latency=winstim(kk,:);
        cfg.covariancewindow=winstim(kk,:);
        timefreqstim{kk,1}=ft_timelockanalysis(cfg,freqdata); % can use to filter also, over whole trial first
        cfg.trials=find(timefreqstim{kk}.usetrial);
        cfg.latency=winstim(1,:);
        cfg.covariancewindow=winstim(1,:);
        timefreqstim{1,kk}=ft_timelockanalysis(cfg,freqdata); % can use to filter also, over whole trial first
        cfg=[];
        cfg.vol=ft_headshape;
        cfg.grid=grid;
        cfg.method='lcmv';
        cfg.keepfilter='yes';
        cfg.fixedori='no';
        source{kk}=ft_sourceanalysis(cfg,timefreqstim{1,kk},timefreqstim{kk,1});
        cfg.grid.filter=source{kk}.filter;
        weights.filter=cfg.grid.filter;
        weights.usetrial=find(timefreqstim{kk,1}.usetrial);
        save(['weights_stimlock_' num2str(freq(jj,1)) 'to' num2str(freq(jj,2)) 'Hz_' num2str(1000*winstim(kk,1)) 'to' num2str(1000*winstim(kk,2)) 'ms_LCMV.mat'],'weights');
        clear weights;
        cfg.keeppow='yes'; % jz added this
        cfg.keepfilter='no';
        cfg.lcmv.keepmom='no';
        cfg.lcmv.powmethod='lambda1'; % or trace?
        sourcecon{kk}=ft_sourceanalysis(cfg,timefreqstim{1,kk});
        sourceact{kk}=ft_sourceanalysis(cfg,timefreqstim{kk,1});
        beamsa=nut_ft2beam(sourceact{kk},mriav2r,'dataset02.hdm',flags);
        beamsc=nut_ft2beam(sourcecon{kk},mriav2r,'dataset02.hdm',flags);
        beam=beamsa;
        beam.s{2}=beamsc.s{1};
        save(['s_beamtf_stimlock_' num2str(freq(jj,1)) 'to' num2str(freq(jj,2)) 'Hz_' num2str(1000*winstim(kk,1)) 'to' num2str(1000*winstim(kk,2)) 'ms_LCMV.mat'],'-struct','beam');
    end
    clear timefreqstim
    clear source*
    clear freqdata
end
nut_tfbf2timef('stimlock','LCMV');


%% response-locked analysis
for jj=1:5
    load(['freqdata' num2str(jj) '.mat'])
    flags.avetime=1;
    fromresp=[-[2.1:-.1:-0.4]; -[1.7:-.1:-.8]]';
    conind=size(fromresp,1);
    for kk=[conind size(fromresp,1)-2:-1:1]
        tfresp{kk}=[];
        for ll=1:length(freqdata.time)
            winresp=[freqdata.time{ll}(end)-1-[2.1:-.1:-0.4]; freqdata.time{ll}(end)-1-[1.7:-.1:-.8]]';
            cfg=[];
            cfg.covariance='yes';
            cfg.permutation='no';
            cfg.keeptrials='no';
            cfg.latency=winresp(kk,:);
            cfg.vartrllength=1;
            cfg.covariancewindow=winresp(kk,:);
            cfg.trials=ll;
            cfg.keeptrials='no';
            cfg.covariance='yes';
            timefreqresp{ll}=ft_timelockanalysis(cfg,freqdata);
            if any(timefreqresp{ll}.usetrial)
                if isempty(tfresp{kk})
                    tfresp{kk}.fsample=timefreqresp{ll}.fsample;
                    tfresp{kk}.label=timefreqresp{ll}.label;
                    tfresp{kk}.dimord=timefreqresp{ll}.dimord;
                    tfresp{kk}.grad=timefreqresp{ll}.grad;
                    tfresp{kk}.time=timefreqresp{ll}.time;
                    
                end
                tfresp{kk}.cov(ll,:,:)=timefreqresp{ll}.cov;
                tfresp{kk}.trialinfo(ll,:)=timefreqresp{ll}.trialinfo;
            else
                tfresp{kk}.cov(ll,:,:)=NaN;
                tfresp{kk}.trialinfo(ll,:)=NaN;
            end
            tfresp{kk}.usetrial(ll)=logical(timefreqresp{ll}.usetrial);
        end
        timefreqresp=[];
    end
    
    for kk=[size(fromresp,1)-2:-1:1]
        cfg=[];
        cfg.vol=ft_headshape;
        cfg.grid=grid;
        cfg.method='lcmv';
        cfg.keepfilter='yes';
        cfg.fixedori='no';
        usetemp=find(~isnan(tfresp{kk}.cov(:,100,100)));
        conuse=tfresp{conind};
        conuse.cov=squeeze(mean(tfresp{conind}.cov(usetemp,:,:),1));
        actuse=tfresp{kk};
        actuse.cov=squeeze(mean(tfresp{kk}.cov(usetemp,:,:),1));
        conuse.avg=zeros(183,2);
        actuse.avg=zeros(183,2);
        sourcer{kk}=ft_sourceanalysis(cfg,conuse,actuse);
        cfg.grid.filter=sourcer{kk}.filter;
        weights.filter=cfg.grid.filter;
        weights.trluse=find(tfresp{kk}.usetrial);
        save(['weights_resplock_' num2str(freq(jj,1)) 'to' num2str(freq(jj,2)) 'Hz_' num2str(1000*fromresp(kk,1)) 'to' num2str(1000*fromresp(kk,2)) 'ms_LCMV.mat'],'weights');
        clear weights
        cfg.keeppow='yes'; % jz added this
        cfg.keepfilter='no';
        cfg.lcmv.keepmom='no';
        cfg.lcmv.powmethod='lambda1'; % or trace?
        sourcercon=ft_sourceanalysis(cfg,conuse);
        sourceract=ft_sourceanalysis(cfg,actuse);
        beamra=nut_ft2beam(sourceract,mriav2r,'dataset02.hdm',flags);
        beamrc=nut_ft2beam(sourcercon,mriav2r,'dataset02.hdm',flags);
        beam=beamra;
        beam.s{2}=beamrc.s{1};
        save(['s_beamtf_resplock_' num2str(freq(jj,1)) 'to' num2str(freq(jj,2)) 'Hz_' num2str(1000*fromresp(kk,1)) 'to' num2str(1000*fromresp(kk,2)) 'ms_LCMV.mat'],'-struct','beam')
    end
    clear tfresp
    clear source*
    clear freqdata
end
nut_tfbf2timef('resplock','LCMV');


