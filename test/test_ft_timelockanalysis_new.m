function test_ft_timelockanalysis_new(datainfo,writeflag)

% MEM 1500mb
% WALLTIME 00:10:00

% ft_timelockanalysis_new ft_timelockanalysis ref_datasets

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

% this is a function for testing ft_timelockanalysis_new, which is not official yet

% the optional writeflag determines whether the output should be saved to
% disk

%% 
% This function is testing a new ft_timelockanalysis_new which is something
% Johanna is working on and not in SVN yet.
return;

%%

if nargin<2
  writeflag = 0;
end
if nargin<1
  datainfo = ref_datasets;
end
% for k = 1:numel(datainfo)
for k = 1:10
  datanew = timelockanalysis10trials(datainfo(k), writeflag);
  
  fname = fullfile(datainfo(k).origdir,'latest/timelock',datainfo(k).type,'timelock_',datainfo(k).datatype);
  tmp = load(fname);
  if isfield(tmp, 'data')
    data = tmp.data;
  elseif isfield(tmp, 'datanew')
    data = tmp.datanew;
  else isfield(tmp, 'timelock')
    data = tmp.timelock;
  end
  
  datanew = rmfield(datanew, 'cfg'); % these are per construction different if writeflag = 0;
  data    = rmfield(data,    'cfg');
  assert(isequaln(data, datanew));
end

test_cfg_options;

function [tlck] = timelockanalysis10trials(dataset, writeflag)

cfg        = [];
cfg.inputfile  = fullfile(dataset.origdir,'latest/raw',dataset.type,['preproc_',dataset.datatype]);
if writeflag
  cfg.outputfile = fullfile(dataset.origdir,'latest/timelock',dataset.type,'timelock_',dataset.datatype);
end
tlck = ft_timelockanalysis(cfg);
tlck1 = ft_timelockanalysis_new(cfg);
return


function test_cfg_options
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/eventrelatedaveraging/dataFC_LP.mat'));
data = dataFC_LP;
clear dataFC_LP

data.time{2}=data.time{2}+.5; % purposely add some jitter to time window
data.time{3}=data.time{3}-.5;

cfg=[];
try
    tlock=ft_timelockanalysis_new(cfg,data);
catch me
%     if ~strcmp(me.message,'the option "output" was not specified or was empty');
        error(me.message)
%     end
end
cfg=[];
cfg.output='rubbish';
try
    tlock=ft_timelockanalysis_new(cfg,data);
catch me
    if ~strcmp(me.message,'the value of cfg.output is not set correctly');
        error(me.message)
    end
end


% no latency or covlatency given, use defaults
cfg=[];
cfg.output='avg';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);
cfg=[];
cfg.output='cov';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.keeptrials='no';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg=[];
cfg.output='avgandcov';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.keeptrials='no';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
% cfg.covariance='yes';
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);
cfg=[];
cfg.output='avg';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.keeptrials='yes';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);
cfg=[];
cfg.output='cov';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.keeptrials='yes';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg=[];
cfg.output='avgandcov';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.keeptrials='yes';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
% cfg.covariance='yes';
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);

% options with .latency and .covlatency specified
cfg=[];
cfg.output='avg';
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.keeptrials='no';
cfg.latency=[min(data.time{1}) max(data.time{1})];
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);
cfg.latency='maxperlength';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='minperlength';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='prestim';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='poststim';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='avg';
cfg.keeptrials='yes';
cfg.latency=[min(data.time{1}) max(data.time{1})];
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='maxperlength';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='minperlength';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='prestim';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
cfg.latency='poststim';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='avgandcov';
cfg.keeptrials='no';
cfg.latency=[min(data.time{1}) max(data.time{1})];
cfg.covlatency=[min(data.time{1}) max(data.time{1})];
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
% cfg.covariance='yes';
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);
cfg.latency='maxperlength';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end


cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='avgandcov';
cfg.keeptrials='yes';
cfg.latency=[min(data.time{1}) max(data.time{1})];
cfg.covlatency=[min(data.time{1}) max(data.time{1})];
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
% cfg.covariance='yes';
% cfg.vartrllength=1;
% tlocko=ft_timelockanalysis(cfg,data);

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.keeptrials='no';
cfg.latency=[min(data.time{1}) max(data.time{1})];
cfg.covlatency=[min(data.time{1}) max(data.time{1})];
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.keeptrials='yes';
cfg.latency=[min(data.time{1}) max(data.time{1})];
cfg.covlatency=[min(data.time{1}) max(data.time{1})];
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

% check toi options
cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.latency=[min(data.time{1}) max(data.time{1})];
cfg.covlatency='minperlength';
cfg.toi=[-0.5 0.7];
cfg.timwin=1;
cfg.equatenumtrials='yes';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.covlatency='minperlength';
cfg.toi=[-.8:.3:0.1];
cfg.timwin=1;
cfg.equatenumtrials='no';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.covlatency='minperlength';
cfg.toi=[-.8:.3:0.1];
cfg.timwin=1;
cfg.equatenumtrials='no';
cfg.keeptrials='yes';
try
    tlock=ft_timelockanalysis_new(cfg,data);
    if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end
catch me
    if ~strcmp(me.message,'sorry, if keeping trials and computing cov, cfg.equatenumtrials should be yes')
        error(me.message)
    end
end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.covlatency='minperlength';
cfg.toi=[-.8:.3:0.1];
cfg.timwin=1;
cfg.equatenumtrials='yes';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

cfg=[];
cfg.feedback='none';
cfg.preproc.feedback='textbar';
cfg.output='cov';
cfg.covlatency='minperlength';
cfg.toi=[-.8:.3:0.1];
cfg.timwin=1;
cfg.equatenumtrials='yes';
cfg.keeptrials='yes';
tlock=ft_timelockanalysis_new(cfg,data);
if ~strmatch(ft_datatype(tlock),'timelock'),error('datatype');end;
if ~strcmp(ft_datatype(ft_checkdata(tlock,'datatype','raw')),'raw') || ~strcmp(ft_datatype(ft_checkdata(tlock)),'timelock'),error('checkdata');end

%% ignore for svn testing, but used to test output of new function in further functions

% test ft_sourceanalysis versus MNE event related tutorial
if 0
    
    cfg=[];
    cfg.output='avg';
    cfg.feedback='none';
    cfg.preproc.feedback='textbar';
    tlock=ft_timelockanalysis_new(cfg,data);
    cfg=[];
    cfg.method='lcmv';
    cfg.hdmfile=dccnpath('/home/common/matlab/fieldtrip/data/Subject01.hdm');
    cfg.grad=data.grad;
    source=ft_sourceanalysis(cfg,tlock);
    
    
    cfg = [];
    cfg.covariance = 'yes';
    cfg.vartrllength=1;
    cfg.covariancewindow = [-inf 0]; %it will calculate the covariance matrix
    % on the timepoints that are
    % before the zero-time point in the trials
    tlckFC = ft_timelockanalysis(cfg, data);
%     tlckFIC = ft_timelockanalysis(cfg, dataFIC_LP);
    save tlck tlckFC tlckFIC;
    cfg=[];
    cfg.method='lcmv';
    cfg.hdmfile=dccnpath('/home/common/matlab/fieldtrip/data/Subject01.hdm');
    cfg.grad=data.grad;
    sourceFC=ft_sourceanalysis(cfg,tlckFC);
end

% spin off of beamformer tutorial but in time-domain, results won't match exactly
if 0
    cfg=[];
    cfg.output='cov';
    cfg.feedback='none';
    cfg.preproc.feedback='textbar';
    cfg.covlatency=[0.8 1.3];
    cfg.preproc.bpfilter='yes';
    cfg.preproc.bpfreq=[16 20];
    tlock=ft_timelockanalysis_new(cfg,data);
    load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'))
    cfg=[];
    vol=ft_prepare_singleshell(cfg,segmentedmri);
    cfg=[];
    cfg.vol=vol;
    cfg.reducerank = 2;
    cfg.grad=tlock.grad;
    cfg.grid.resolution=1;
    cfg.channel={'MEG','-MLP31','-MLO12'};
    grid=ft_prepare_leadfield(cfg);
    cfg=[];
    cfg.method='lcmv';
    cfg.projectnoise='yes';
    cfg.grid=grid;
    cfg.vol=vol;
    source=ft_sourceanalysis(cfg,tlock);
    
    mri=ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));
    sourcediff=source;
    % sourcediff.avg.pow=(source.avg.pow-source.avg.noise)./source.avg.noise;
    sourcediff.avg.pow=(source.avg.pow)./source.avg.noise;
    cfg=[];
    cfg.downsample=2;
    sourcediffint=ft_sourceinterpolate(cfg,sourcediff,mri);
    cfg=[];
    cfg.method='slice';
    cfg.funparameter='avg.pow';
    cfg.maskparameter=cfg.funparameter;
    cfg.funcolorlim=[5 6.2];
    cfg.opacitylim=[5 6.2];
    cfg.opacitymap='rampup';
    figure;ft_sourceplot(cfg,sourcediffint); %ok
    
    cfg=[];
    cfg.downsample=2;
    sourceint=ft_sourceinterpolate(cfg,source,mri);
    cfg=[];
    cfg.method='slice';
    cfg.funparameter='avg.pow';
    % cfg.maskparameter=cfg.funparameter;
    % cfg.funcolorlim=[5 6.2];
    % cfg.opacitylim=[5 6.2];
    % cfg.opacitymap='rampup';
    figure;ft_sourceplot(cfg,sourceint); %ok
    
    
end

% test ft_timelockstatistics
if 0
    
cfg=[];
cfg.output='cov';
cfg.feedback='none';
cfg.keeptrials='yes';
cfg.preproc.feedback='textbar';
tlock=ft_timelockanalysis_new(cfg,data);
    
tlock1 = ft_selectdata(tlock, 'rpt',1:36)
tlock2 = ft_selectdata(tlock, 'rpt',37:72)

cfg=[];
cfg.method='analytic';
   cfg.statistic='ft_statfun_indepsamplesT' ;
stat=ft_timelockstatistics(cfg,tlock1,tlock2);
end

