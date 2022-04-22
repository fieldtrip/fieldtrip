function test_example_common_filters_in_beamforming

% MEM 4gb
% WALLTIME 00:10:00

%
%% Common filters in beamforming
%
% When you want to reconstruct and compare the sources of two (or more) conditions using a beamformer approach, you can either compute separate spatial filters for each condition, or use a common filter based on the combined conditions.
%
%% # What are the advantages of common filters?
%
% The main advantage of using common filters is that more data is used for constructing the spatial filters. The data of the conditions is appended, and the cross-spectral density matrix is based on this combined dataset. Generally, the more data is used in this step, the better the estimate, thus the more reliable the filters.
% Another advantage is that you use the same filters for estimating the sources in both conditions. Differences in source activity can then be ascribed to power differences in conditions, not to differences between the filters.
%
%% # When can you use common filters?
%
% Common filters can be used when you want to compare conditions for which you assume the underlying sources are the same, but active to a different extent.
% It is not a problem to have different amounts of trials in the conditions.
% It is a requirement that the time windows in all conditions are of equal length.
%
%% # How to do this in FieldTrip
%
% In FieldTrip, common filters can be used both with the DICS and the PCC beamformer approach. If you are interested in calculating source reconstructed power, both methods can be used and will lead to similar results. PCC is faster but more memory-demanding, whereas DICS is slower but more memory-friendly. The choice of methods depends on personal preferences, your data (e.g., number of trials, number of tapers) and computer specs.
%
% The general procedure is as follows:
%
% #  calculate the cross-spectral density matrix of the combined conditions
% #  compute the spatial filters
% #  project all single trials through these filters
% #  compute average for each condition (or do statistics)
%
%% # Example code
%
% FIXME Code for reconstruction of single trial data is incomplete
%
% In the following, scripts for both approaches are presented.
% Let's say we have data for two conditions, condition A and B, and we assume that the same sources are active in both, but to a different extent.
%
% We have the preprocessed data for both conditions (_dataA_ and _dataB_) and we precomputed the sourcemodel and the volume conduction model (_sourcemodel_ and _vol_). For more information on how to do this, have a look at the [beamformer tutorial](/tutorial/beamformer) and the FieldTrip functions **[ft_prepare_leadfield](https://github.com/fieldtrip/fieldtrip/blob/release/ft_prepare_leadfield.m)** and **[ft_prepare_headmodel](https://github.com/fieldtrip/fieldtrip/blob/release/ft_prepare_headmodel.m)** with method='singleshell'.
%
%% ## PCC
%
% append the two conditions and remember the design %

% JM added for this function: load some data
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataPre.mat')); dataA = dataPre; clear dataPre;
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/dataPost.mat')); dataB = dataPost; clear dataPost;
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/headmodel.mat'));
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/sourcemodel.mat'));


data = ft_appenddata([], dataA, dataB);
design = [ones(1,length(dataA.trial)) ones(1,length(dataB.trial))*2]; % only necessary if you are interested in reconstructing single trial data

% ft_freqanalysis %
cfg=[];
cfg.method      = 'mtmfft';
cfg.output      = 'fourier';  % gives the complex Fourier spectra
cfg.foilim      = [60 60];    %analyze40-80 Hz (60 Hz +/- 20 Hz smoothing)
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 20;
cfg.keeptrials  = 'yes';      % in order to separate the conditions again afterwards, we need to keep the trials. This is not otherwise necessary to compute the common filter
cfg.keeptapers  = 'yes';

freq = ft_freqanalysis(cfg, data);

% compute common spatial filter AND project all trials through it %
cfg=[];
cfg.method      = 'pcc';
cfg.sourcemodel = sourcemodel;       % previously computed sourcemodel
cfg.headmodel   = headmodel;        % previously computed volume conduction model
cfg.frequency   = 60;
cfg.keeptrials  = 'yes';      % keep single trials. Only necessary if you are interested in reconstructing single trial data

source = ft_sourceanalysis(cfg, freq);

% average over tapers, keep single trials %

% This step is only necessary if you need to reconstruct single trial data

cfg=[];
cfg.keeptrials = 'yes';

source = ft_sourcedescriptives(cfg, source); % contains the source estimates for all trials/both conditions

% calculate average for each condition %
A = find(design==1); % find trial numbers belonging to condition A
B = find(design==2); % find trial numbers belonging to condition B

sourceA = source;
sourceA.trial(B) = [];
sourceA.cumtapcnt(B) = [];
sourceA.method = 'rawtrial';
sourceA = ft_sourcedescriptives([], sourceA); % compute average source reconstruction for condition A

sourceB=source;
sourceB.trial(A) = [];
sourceB.cumtapcnt(A) = [];
sourceB.method = 'rawtrial';
sourceB = ft_sourcedescriptives([], sourceB); % compute average source reconstruction for condition B

%% ## DICS
%
% append the two conditions and remember the design %
data = ft_appenddata([], dataA, dataB);
design = [ones(1,length(dataA.trial)) ones(1,length(dataB.trial))*2]; % only necessary if you are interested in reconstructing single trial data

% freqanalysis %
cfg=[];
cfg.method      = 'mtmfft';
cfg.output      = 'powandcsd';  % gives power and cross-spectral density matrices
cfg.foilim      = [60 60];      %analyze40-80 Hz (60 Hz +/- 20 Hz smoothing)
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 20;
cfg.keeptrials  = 'yes';        % in order to separate the conditions again afterwards, we need to keep the trials. This is not otherwise necessary to compute the common filter
cfg.keeptapers  = 'no';

freq = ft_freqanalysis(cfg, data);

% compute common spatial filter %
cfg=[];
cfg.method      = 'dics';
cfg.sourcemodel = sourcemodel;         % previously computed sourcemodel
cfg.headmodel   = headmodel;          % previously computed volume conduction model
cfg.frequency   = 60;
cfg.dics.keepfilter  = 'yes';        % remember the filter

source = ft_sourceanalysis(cfg, freq);

% project all trials through common spatial filter %
cfg=[];
cfg.method      = 'dics';
cfg.sourcemodel = sourcemodel;       % previously computed grid
cfg.headmodel   = headmodel;        % previously computed volume conduction model
cfg.sourcemodel.filter = source.avg.filter; % use the common filter computed in the previous step!
cfg.frequency   = 60;
cfg.rawtrial    = 'yes';      % project each single trial through the filter. Only necessary if you are interested in reconstructing single trial data

source = ft_sourceanalysis(cfg, freq); % contains the source estimates for all trials/both conditions

% calculate average for each condition %

% This step is only necessary if you need to reconstruct single trial data

A = find(design==1); % find trial numbers belonging to condition A
B = find(design==2); % find trial numbers belonging to condition B

sourceA = source;
sourceA.trial(B) = [];
sourceA.cumtapcnt(B) = [];
sourceA.df = length(A);
sourceA = ft_sourcedescriptives([], sourceA); % compute average source reconstruction for condition A

sourceB=source;
sourceB.trial(A) = [];
sourceB.cumtapcnt(A) = [];
sourceB.df = length(B);
sourceB = ft_sourcedescriptives([], sourceB); % compute average source reconstruction for condition B
