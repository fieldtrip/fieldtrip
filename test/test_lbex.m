function test_lbex(datadirs)

% WALLTIME 01:00:00
% MEM 12gb
% DEPENDENCY ft_prepare_leadfield

if nargin==0
  datadirs{1} = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf');
  datadirs{2} = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/networkanalysis');
end

%% read the continuous data and segment into 2 seconds epochs, with 50% overlap
cfg            = [];
cfg.dataset    = fullfile(datadirs{1},'SubjectRest.ds');
cfg.continuous = 'yes';
cfg.channel    = {'MEG'};
data1 = ft_preprocessing(cfg);
cfg.coilaccuracy = 1;
data2 = ft_preprocessing(cfg);

cfg         = [];
cfg.length  = 2;
cfg.overlap = 0.5;
data1       = ft_redefinetrial(cfg, data1);
data2       = ft_redefinetrial(cfg, data2);

cfg        = [];
cfg.demean = 'yes';
cfg.trials = 1:(numel(data1.trial)-6);
data1      = ft_preprocessing(cfg, data1);
data2      = ft_preprocessing(cfg, data2);

data1 = removefields(data1, 'elec');
data2 = removefields(data2, 'elec');

datadir = datadirs{2};
cd(datadir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BYPASS THE INTERACTIVE PART, DECLARE TRIALS 38 91 153 AS BAD
%% make a visual inspection and reject bad trials/sensors
% cfg = [];
% cfg.method  = 'summary';
% cfg.channel = 'MEG';
% cfg.layout  = 'CTF275.lay';
% dataclean = ft_rejectvisual(cfg, data);
%
%% you can check the rejected trial numbers by typing
% trlind = [];
% for i=1:length(dataclean.cfg.artfctdef.summary.artifact)
%  trlind(i) = find(data.sampleinfo(:,1)==dataclean.cfg.artfctdef.summary.artifact(i));
% end
% disp(trlind);

badtrials  = [18 19 21 72 73 74 75 76 93 94 109 110 126 127 128 140 172 173 179 180 181 182 196 197 198 227 228 233 243 244 250 251 265 266 286];
cfg        = [];
cfg.trials = setdiff(1:numel(data1.trial), badtrials);
data1  = ft_selectdata(cfg, data1);
data2  = ft_selectdata(cfg, data2);

%% downsample the data to speed up component analysis
data1.time(1:end) = data1.time(1);
data2.time(1:end) = data2.time(1);

cfg            = [];
cfg.resamplefs = 100;
cfg.detrend    = 'yes';
data1         = ft_resampledata(cfg, data1);
data2         = ft_resampledata(cfg, data2);

%% use ICA in order to identify cardiac and blink components
%cfg                 = [];
%cfg.method          = 'runica';
%cfg.runica.maxsteps = 50;
%cfg.randomseed      = 0;
%comp                = ft_componentanalysis(cfg, datads);
load(fullfile(datadir,'comp.mat'));

%% visualize components

% these were the indices of the bad comp**[[reference:ft_definetrial|ft_definetrial]]** and onents that were identified
% they may be different if you re-run the ICA decomposition
badcomp = [2 3 7 16];

cfg           = [];
cfg.component = badcomp;
data1      = ft_rejectcomponent(cfg, comp, data1);
data2      = ft_rejectcomponent(cfg, comp, data2);

%% load the required geometrical information
load(fullfile(datadir, 'hdm.mat'));

data1.grad = ft_convert_units(data1.grad, 'm');
hdm = ft_convert_units(hdm, 'm');

%% compute sourcemodels
cfg            = [];
cfg.method     = 'basedonresolution';
cfg.resolution = 0.004;
cfg.headmodel  = hdm;
sourcemodel4   = ft_prepare_sourcemodel(cfg);
%cfg.resolution = 0.008;
%sourcemodel8   = ft_prepare_sourcemodel(cfg);

sel = false(sourcemodel4.dim);
sel(1:2:end,1:2:end,1:2:end) = true;
sourcemodel8 = sourcemodel4;
sourcemodel8.pos = sourcemodel8.pos(sel(:),:);
sourcemodel8.inside = sourcemodel8.inside(sel(:));
sourcemodel8.dim = ceil(sourcemodel4.dim/2);

%% compute the leadfield
cfg             = [];
cfg.sourcemodel = sourcemodel8;
cfg.headmodel   = hdm;
cfg.channel     = {'MEG'};
cfg.singleshell.batchsize = 1500;
lf81            = ft_prepare_leadfield(cfg, data1);
lf82            = ft_prepare_leadfield(cfg, data2);

cfg.lbex = 0.02;
lf81lbex = ft_prepare_leadfield(cfg, data1);
lf82lbex = ft_prepare_leadfield(cfg, data2);

%crank up the order
cfg = rmfield(cfg, 'lbex');
cfg.headmodel.order = 20;
lf8120          = ft_prepare_leadfield(cfg, data1);
lf8220          = ft_prepare_leadfield(cfg, data2);

cfg.lbex = 0.025;
lf8120lbex = ft_prepare_leadfield(cfg, data1);
lf8220lbex = ft_prepare_leadfield(cfg, data2);

cfg.sourcemodel = sourcemodel4;
cfg.keep   = sel(:);
lf4120lbex = ft_prepare_leadfield(cfg, data1);
lf4220lbex = ft_prepare_leadfield(cfg, data2);

prune = {'leadfield' 'subspace'};% 'inside' 'pos'};
for k = 1:numel(prune)
  lf4120lbex.(prune{k}) = lf4120lbex.(prune{k})(sel(:));
  lf4220lbex.(prune{k}) = lf4220lbex.(prune{k})(sel(:));
end
prune = {'inside' 'pos'};
for k = 1:numel(prune)
  lf4120lbex.(prune{k}) = lf4120lbex.(prune{k})(sel(:),:);
  lf4220lbex.(prune{k}) = lf4220lbex.(prune{k})(sel(:),:);
end


%% compute sensor level Fourier spectra
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;
freq           = ft_freqanalysis(cfg, data1);


%% compute the actual source reconstruction
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'dics';
cfg.headmodel         = hdm;
cfg.keeptrials        = 'yes';
cfg.dics.lambda       = '20%';
cfg.dics.projectnoise = 'yes';
cfg.dics.fixedori     = 'yes';
cfg.dics.realfilter   = 'yes';
cfg.dics.keepfilter   = 'yes';
cfg.sourcemodel       = lf81;
source81              = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf82;
source82              = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf8120;
source8120            = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf8220;
source8220            = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf8120lbex;
source8120lbex        = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf8220lbex;
source8220lbex        = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf4120lbex;
source4120lbex        = ft_sourceanalysis(cfg, freq);
cfg.sourcemodel       = lf4220lbex;
source4220lbex        = ft_sourceanalysis(cfg, freq);

freqC = ft_checkdata(freq, 'cmbrepresentation', 'fullfast');
C     = freqC.crsspctrm;

F81 = cat(1,source81.avg.filter{:});
C81 = F81*C*F81';
C81 = C81./sqrt(diag(C81)*diag(C81)');

F82 = cat(1,source82.avg.filter{:});
C82 = F82*C*F82';
C82 = C82./sqrt(diag(C82)*diag(C82)');

F8120 = cat(1,source8120.avg.filter{:});
C8120 = F8120*C*F8120';
C8120 = C8120./sqrt(diag(C8120)*diag(C8120)');

F8220 = cat(1,source8220.avg.filter{:});
C8220 = F8220*C*F8220';
C8220 = C8220./sqrt(diag(C8220)*diag(C8220)');

F8120lbex = cat(1,source8120lbex.avg.filter{:});
C8120lbex = F8120lbex*C*F8120lbex';
C8120lbex = C8120lbex./sqrt(diag(C8120lbex)*diag(C8120lbex)');

F8220lbex = cat(1,source8220lbex.avg.filter{:});
C8220lbex = F8220lbex*C*F8220lbex';
C8220lbex = C8220lbex./sqrt(diag(C8220lbex)*diag(C8220lbex)');

F4120lbex = cat(1,source4120lbex.avg.filter{:});
C4120lbex = F4120lbex*C*F4120lbex';
C4120lbex = C4120lbex./sqrt(diag(C4120lbex)*diag(C4120lbex)');

F4220lbex = cat(1,source4220lbex.avg.filter{:});
C4220lbex = F4220lbex*C*F4220lbex';
C4220lbex = C4220lbex./sqrt(diag(C4220lbex)*diag(C4220lbex)');

% conclusion: lbex leads to overall less spurious coherence
figure; histogram(abs(C8120(:)),20)
hold on; histogram(abs(C4120lbex(:)),20)

