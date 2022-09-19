function test_dss_ecg

% MEM 8gb
% WALLTIME 00:10:00
% DEPENDENCY ft_artifact_zvalue ft_componentanalysis

% define 'trl' of example dataset and read in the meg data
cfg         = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/ArtifactRemoval.ds');
cfg.trialdef.eventtype = 'trial';
cfg     = ft_definetrial(cfg);
cfg.trl = cfg.trl(1:end-1,:); % remove the last one, because it clips.

cfg.channel = 'MEG';
cfg.demean  = 'yes';
meg         = ft_preprocessing(cfg);


cfg = removefields(cfg, {'channel', 'demean'});
cfg.artfctdef.zvalue.channel         = 'EEG058';
cfg.artfctdef.zvalue.cutoff          = 2;
%cfg.artfctdef.zvalue.interactive     = 'yes';
cfg.artfctdef.zvalue.bpfilter        = 'yes';
cfg.artfctdef.zvalue.bpfreq          = [5 30];
cfg.artfctdef.zvalue.hilbert         = 'yes';
cfg.artfctdef.zvalue.artfctpeak      = 'yes';
cfg.artfctdef.zvalue.artfctpeakrange = [-.25 .5]; % save out 250ms prior and 500ms post ECG peak
cfg = ft_artifact_zvalue(cfg);


% paramscell.tr  = cfg.artfctdef.zvalue.peaks_indx;
% paramscell.pre = 0.25*meg.fsample;
% paramscell.pst = 0.50*meg.fsample;
% paramscell.demean = true;
paramscell.artifact = cfg.artfctdef.zvalue.artifact;

cfg                   = [];
cfg.method            = 'dss';
cfg.dss.denf.function = 'denoise_avg2';
cfg.dss.denf.params   = paramscell;
cfg.dss.wdim          = 75;
cfg.numcomponent      = 4;
cfg.channel           = 'MEG';
cfg.cellmode          = 'yes';
comp                  = ft_componentanalysis(cfg, meg);

%cfg = [];
%cfg.layout = 'CTF275_helmet.mat';
%ft_databrowser(cfg, comp);

cfg = [];
cfg.component = [1 2];
meg_clean     = ft_rejectcomponent(cfg, comp, meg);
