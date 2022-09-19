function test_pull1377

% MEM 12gb
% WALLTIME 1:30:00
% DEPENDENCY ft_prepare_sourcemodel headsurface ft_prepare_leadfield ft_freqanalysis ft_sourceanalysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function creates a set of input-structures to be used for testing source
% analysis steps prior to inverse solution. The structure of this script is more or less
% 1. get the data (MEG and EEG)
% 2. preprocess the data (MEG and EEG)
%   a. timelock with covariance for lcmv and mne
%   b. timelock without keeptrials for dipolefitting
%   c. freq for dics and pcc
% 3. create leadfield
%   a. internally computed leadfield (prepare volume conductor first and leadfield afterwards (MEG localspheres, EEG single sphere))
%   b. load externally computed leadfield and convert to appropriate format
% 4. perform source analysis (MEG and EEG): LCMV, MNE, DICS, PCC, dipolefitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In general, variables in one function workspace are not available to other
% functions. However, nested functions can access and modify variables in the
% workspaces of the functions that contain them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% prevent errors from cfg.mne.keepleadfield, etc
global ft_default
ft_default.checkconfig = 'loose';


datameg = []; % touch these variables to make them shared between the main function and the nexted function
dataeeg = [];

loading_raw_data

headmodel_megsphere = [];
headmodel_eegsphere = [];
sourcemodel_megsphere = [];
sourcemodel_eegsphere = [];
sourcemodel_precomputed_meg = [];
sourcemodel_precomputed_eeg = [];
headmodel_simbio = [];
sourcemodel_simbio = [];
headmodel_openmeeg = [];
sourcemodel_openmeeg = [];

create_leadfields

MEG_tlck = [];
MEG_tlck_df = [];
MEG_freq = [];
EEG_tlck = [];
EEG_tlck_df = [];
EEG_freq = [];

preprocessing

sourcemodel_meg = sourcemodel_megsphere;
headmodel_meg = headmodel_megsphere;
sourcemodel_eeg = sourcemodel_eegsphere;
headmodel_eeg = headmodel_eegsphere;
perform_source_analysis

sourcemodel_meg = sourcemodel_precomputed_meg;
headmodel_meg = [];
sourcemodel_eeg = sourcemodel_precomputed_eeg;
headmodel_eeg = [];
perform_source_analysis

sourcemodel_meg = [];
headmodel_meg = [];
sourcemodel_eeg = sourcemodel_simbio;
headmodel_eeg = headmodel_simbio;
perform_source_analysis

sourcemodel_meg = [];
headmodel_meg = [];
sourcemodel_eeg = sourcemodel_openmeeg;
headmodel_eeg = headmodel_openmeeg;
perform_source_analysis

% !! create worse-case scenario, whereby order and nr of chans don't match across inputs (grad/elec and leadfields)
% Step 1: remove 2-3 random chans from both MEG and EEG raw data
% Step 2: shuffle the order of the channels (note that this breaks the mapping between channels and data)
cfg = [];
cfg.channel = randperm(length(datameg.label)-2);
datameg = ft_selectdata(cfg, datameg);
datameg.label = datameg.label(randperm(length(datameg.label)));

cfg.channel = randperm(length(dataeeg.label)-2);
dataeeg = ft_selectdata(cfg, dataeeg);
dataeeg.label = dataeeg.label(randperm(length(dataeeg.label)));

preprocessing

sourcemodel_meg = sourcemodel_megsphere;
headmodel_meg = headmodel_megsphere;
sourcemodel_eeg = sourcemodel_eegsphere;
headmodel_eeg = headmodel_eegsphere;
perform_source_analysis

sourcemodel_meg = sourcemodel_precomputed_meg;
headmodel_meg = [];
sourcemodel_eeg = sourcemodel_precomputed_eeg;
headmodel_eeg = [];
perform_source_analysis

sourcemodel_meg = [];
headmodel_meg = [];
sourcemodel_eeg = sourcemodel_simbio;
headmodel_eeg = headmodel_simbio;
perform_source_analysis

sourcemodel_meg = [];
headmodel_meg = [];
sourcemodel_eeg = sourcemodel_openmeeg;
headmodel_eeg = headmodel_openmeeg;
perform_source_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. get the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function loading_raw_data
    
    clc
    disp('1: Get the raw data')
    
    % for MEG data + sensor info
    
    load(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/preproc_ctf151.mat'), 'data');
    datameg = data;
    clear data
    
    % for EEG data + channel info
    
    cfg = [];
    cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/preprocessing_erp/s04.eeg');
    hdr = ft_read_header(cfg.dataset);
    event = ft_read_event(cfg.dataset);
    EVsample = [event.sample]';
    EVvalue = {event.value}';
    Word = find(strcmp('S141', EVvalue)==1);
    
    begsample = EVsample(Word) - 10;
    endsample = EVsample(Word) + 500;
    offset = -100*ones(size(endsample));
    cfg.trl = [begsample endsample offset];
    
    dataeeg = ft_preprocessing(cfg);
    
    % Fix EEG channel names, but note that the actual recording was not done with a Biosemi system
    cfg = [];
    cfg.layout = 'biosemi64.lay';
    lay = ft_prepare_layout(cfg);
    dataeeg.label = lay.label(1:end-2);
    clear lay
    
    dataeeg.elec = ft_read_sens(dccnpath('/home/common/matlab/fieldtrip/template/electrode/standard_1020.elc'), 'senstype', 'eeg');
    
  end % loding raw data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. preprocess the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function preprocessing
    
    clc
    disp('2: preprocess the data')
    
    % create timelock structure with covariance for lcmv and mne
    % and timeloch without keeptrials for dipolefitting
    
    % for MEG
    cfg = [];
    cfg.covariance = 'yes';
    cfg.keeptrials = 'yes';
    cfg.channel = 'MEG';
    MEG_tlck = ft_timelockanalysis(cfg, datameg);
    
    cfg = [];
    cfg.covariance = 'yes';
    % cfg.keeptrials = 'yes'; %if this is not commented, the .avg field necessary in dipolefitting is missing
    cfg.channel = 'MEG';
    MEG_tlck_df = ft_timelockanalysis(cfg, datameg);
    
    % for EEG
    cfg = [];
    cfg.covariance = 'yes';
    cfg.keeptrials = 'yes';
    EEG_tlck = ft_timelockanalysis(cfg, dataeeg);
    
    % for EEG (dimension problems in dipolefitting!!)
    cfg = [];
    cfg.covariance = 'yes';
    % cfg.keeptrials = 'yes';
    EEG_tlck_df = ft_timelockanalysis(cfg, dataeeg);
    
    % create freq structure for dics beamforming and pcc beamforming
    % for MEG
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'fourier';
    cfg.tapsmofrq = 4;
    cfg.foilim = [0 20];
    cfg.channel = 'MEG';
    MEG_freq = ft_freqanalysis(cfg, datameg);
    
    % for EEG
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'fourier';
    cfg.tapsmofrq = 4;
    cfg.foilim = [0 20];
    EEG_freq = ft_freqanalysis(cfg, dataeeg);
    
  end % function preprocessing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. create leadfield
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function create_leadfields
    
    clc
    disp('3: create leadfield')
    
    %% 3.a internally
    
    % get volume conductor model
    % for MEG, localsphere
    headmodel_megsphere = ft_read_headmodel(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/vol/Subject01vol_localspheres.mat'));
    
    % for EEG, singlesphere
    headmodel_eegsphere.o = [0 0 2];
    headmodel_eegsphere.r = 8;
    headmodel_eegsphere.c = 1;
    headmodel_eegsphere.unit = 'cm';
    
    % create leadfield (internally to ft)
    % for MEG
    cfg = [];
    cfg.grad = datameg.grad;
    cfg.headmodel = headmodel_megsphere;
    cfg.channel = 'MEG';
    cfg.resolution = 1.5;
    sourcemodel_megsphere = ft_prepare_leadfield(cfg);
    
    cfg = [];
    cfg.elec = dataeeg.elec;
    cfg.headmodel = headmodel_eegsphere;
    cfg.channel = 'EEG';
    cfg.resolution = 1.5;
    sourcemodel_eegsphere = ft_prepare_leadfield(cfg); % inconsistent with dataeeg.labels
    
    %% 3.b externally (mimic externally created leadfields)
    clc;
    disp('3b: Externally computed leadfields for EEG and MEG')
    
    Ndip = 1000;
    
    % for MEG
    % here we are assuming that the externally computed leadfield has the structure Ndip x Nsens
    leadf = randn(3*Ndip, length(datameg.grad.label));
    
    sourcemodel_precomputed_meg = [];
    sourcemodel_precomputed_meg.leadfield = cell(Ndip,1);
    sourcemodel_precomputed_meg.unit = 'cm';
    sourcemodel_precomputed_meg.leadfielddimord = '{pos}_chan_ori';
    sourcemodel_precomputed_meg.label = datameg.grad.label;
    
    for d =1:Ndip
      sourcemodel_precomputed_meg.leadfield{d} = [leadf(d,:); leadf(d+Ndip,:); leadf(d+2*Ndip,:)]';
    end
    
    sourcemodel_precomputed_meg.inside = ones(size(sourcemodel_precomputed_meg.leadfield));
    sourcemodel_precomputed_meg.pos = randn(length(sourcemodel_precomputed_meg.leadfield),3);
    
    % for EEG
    % here we are assuming that the externally computed leadfield has the structure Ndip x Nsens
    leadf = randn(3*Ndip, length(dataeeg.elec.label));
    
    sourcemodel_precomputed_eeg = [];
    sourcemodel_precomputed_eeg.leadfield = cell(Ndip,1);
    sourcemodel_precomputed_eeg.unit = 'cm';
    sourcemodel_precomputed_eeg.leadfielddimord = '{pos}_chan_ori';
    sourcemodel_precomputed_eeg.label = dataeeg.elec.label;
    
    for d =1:Ndip
      sourcemodel_precomputed_eeg.leadfield{d,1} = [leadf(d,:); leadf(d+Ndip,:); leadf(d+2*Ndip,:)]';
    end
    
    sourcemodel_precomputed_eeg.inside = ones(size(sourcemodel_precomputed_eeg.leadfield));
    sourcemodel_precomputed_eeg.pos = randn(length(sourcemodel_precomputed_eeg.leadfield),3);
    
    %% 3.c simbio leadfield for EEG
    
    segm = [];
    segm.brain = logical(floor(3*rand(10,10,10)));
    segm.dim = size(segm.brain);
    segm.unit = 'mm';
    segm.coordsys = 'ctf';
    segm.transform = eye(4);
    
    cfg = [];
    cfg.method = 'hexahedral';
    cfg.resolution = 1; % this is in mm
    mesh_vol = ft_prepare_mesh(cfg, segm);
    
    ft_plot_mesh(mesh_vol, 'surfaceonly', false);
    view(120, 30)
    
    cfg = [];
    cfg.method ='simbio';
    cfg.conductivity = 0.33; % order follows mesh.tissuelabel
    headmodel_simbio = ft_prepare_headmodel(cfg, mesh_vol);
    
    cfg = [];
    cfg.elec = dataeeg.elec;
    cfg.headmodel = headmodel_simbio;
    cfg.channel = 'EEG';
    cfg.resolution = 1.5;
    sourcemodel_simbio = ft_prepare_leadfield(cfg); % inconsistent with dataeeg.labels
    
    %% 3.c OpenMEEG leadfield for EEG
    
    load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/SubjectSEF_vol.mat'), 'vol')
    bnd_openmeeg = vol;
    
    cfg = [];
    cfg.conductivity = 0.33;
    cfg.method = 'openmeeg';
    headmodel_openmeeg = ft_prepare_headmodel(cfg, bnd_openmeeg);
    
    cfg = [];
    cfg.elec = dataeeg.elec;
    cfg.headmodel = headmodel_openmeeg;
    cfg.channel = 'EEG';
    cfg.resolution = 1.5;
    sourcemodel_openmeeg = ft_prepare_leadfield(cfg); % inconsistent with dataeeg.labels
    
  end % function create_leadfields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. perform source analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function perform_source_analysis
    
    clc;
    disp('4: inverse solution, ft_sourceanalysis and ft_dipolefitting')
    
    %% for MEG
    if ~isempty(sourcemodel_meg) || ~isempty(headmodel_meg)
      
      % do LCMV beamforming
      cfg = [];
      cfg.method = 'lcmv';
      cfg.keepleadfield = 'yes';
      cfg.lcmv.keepfilter = 'yes';
      cfg.lcmv.keepcov = 'yes';
      cfg.lcmv.lambda = '5%';
      cfg.sourcemodel = sourcemodel_meg;
      cfg.headmodel = headmodel_meg;
      sourcelcmv3d1 = ft_sourceanalysis(cfg, MEG_tlck);
      sourcelcmv3d1 = ft_datatype_source(sourcelcmv3d1, 'version', 'latest');
      
      % project through computed filter
      cfg.rawtrial = 'yes';
      cfg.sourcemodel = keepfields(sourcelcmv3d1, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
      ft_sourceanalysis(cfg, MEG_tlck);
      
      % do MNE
      cfg = [];
      cfg.method = 'mne';
      cfg.mne.keepleadfield = 'yes';
      cfg.mne.keepfilter = 'yes';
      cfg.mne.lambda = 1e4;
      cfg.sourcemodel = sourcemodel_meg;
      cfg.headmodel = headmodel_meg;
      sourcemne3d1 = ft_sourceanalysis(cfg, MEG_tlck);
      sourcemne3d1 = ft_datatype_source(sourcemne3d1, 'version', 'latest');
      
      cfg.rawtrial = 'yes';
      cfg.sourcemodel = keepfields(sourcemne3d1, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
      ft_sourceanalysis(cfg, MEG_tlck);
      
      % do DICS
      cfg = [];
      cfg.method = 'dics';
      cfg.dics.keepfilter = 'yes';
      cfg.dics.keepleadfield = 'yes';
      cfg.dics.keepcsd = 'yes';
      cfg.dics.lambda = '5%';
      cfg.frequency = 10;
      cfg.sourcemodel = sourcemodel_meg;
      cfg.headmodel = headmodel_meg;
      sourcedics3d1 = ft_sourceanalysis(cfg, MEG_freq);
      sourcedics3d1 = ft_datatype_source(sourcedics3d1, 'version', 'latest');
      
      cfg.rawtrial = 'yes';
      cfg.sourcemodel = keepfields(sourcedics3d1, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
      ft_sourceanalysis(cfg, MEG_freq);
      
      % do PCC
      cfg = [];
      cfg.method = 'pcc';
      cfg.pcc.keepfilter = 'yes';
      cfg.pcc.keepleadfield = 'yes';
      cfg.pcc.keepcsd = 'yes';
      cfg.pcc.keepmom = 'yes';
      cfg.pcc.lambda = '5%';
      cfg.frequency = 10;
      cfg.sourcemodel = sourcemodel_meg;
      cfg.headmodel = headmodel_meg;
      ft_sourceanalysis(cfg, MEG_freq);
      
      % do dipolefit
      cfg = [];
      cfg.numdipoles = 1;
      % cfg.nonlinear = 'no';
      cfg.sourcemodel = sourcemodel_meg;
      cfg.headmodel = headmodel_meg;
      cfg.latency = 0.025;
      cfg.channel = [1 3:151];
      ft_dipolefitting(cfg, MEG_tlck_df);
    end
    
    %% EEG
    if ~isempty(sourcemodel_eeg) || isempty(headmodel_eeg)
      
      % do LCMV beamforming
      cfg = [];
      cfg.method = 'lcmv';
      cfg.lcmv.keepleadfield = 'yes';
      cfg.lcmv.keepfilter = 'yes';
      cfg.lcmv.keepcov = 'yes';
      cfg.lcmv.lambda = '5%';
      cfg.sourcemodel = sourcemodel_eeg;
      cfg.headmodel = headmodel_eeg;
      sourcelcmv3d1 = ft_sourceanalysis(cfg, EEG_tlck);
      sourcelcmv3d1 = ft_datatype_source(sourcelcmv3d1, 'version', 'latest');
      
      % project through computed filter
      cfg.rawtrial = 'yes';
      cfg.sourcemodel = keepfields(sourcelcmv3d1, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
      ft_sourceanalysis(cfg, EEG_tlck);
      
      % do MNE
      cfg = [];
      cfg.method = 'mne';
      cfg.mne.keepleadfield = 'yes';
      cfg.mne.keepfilter = 'yes';
      cfg.mne.lambda = 1e4;
      cfg.sourcemodel = sourcemodel_eeg;
      cfg.headmodel = headmodel_eeg;
      sourcemne3d1 = ft_sourceanalysis(cfg, EEG_tlck);
      sourcemne3d1 = ft_datatype_source(sourcemne3d1, 'version', 'latest');
      
      cfg.rawtrial = 'yes';
      cfg.sourcemodel = keepfields(sourcemne3d1, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
      ft_sourceanalysis(cfg, EEG_tlck);
      
      % do DICS
      cfg = [];
      cfg.method = 'dics';
      cfg.dics.keepfilter = 'yes';
      cfg.dics.keepleadfield = 'yes';
      cfg.dics.keepcsd = 'yes';
      cfg.dics.lambda = '5%';
      cfg.frequency = 10;
      cfg.sourcemodel = sourcemodel_eeg;
      cfg.headmodel = headmodel_eeg;
      sourcedics3d1 = ft_sourceanalysis(cfg, EEG_freq);
      
      cfg.rawtrial = 'yes';
      cfg.sourcemodel = keepfields(sourcedics3d1, {'pos', 'tri', 'dim', 'inside', 'leadfield', 'leadfielddimord', 'label'});
      ft_sourceanalysis(cfg, EEG_freq);
      
      % do PCC
      cfg = [];
      cfg.method = 'pcc';
      cfg.pcc.keepfilter = 'yes';
      cfg.pcc.keepleadfield = 'yes';
      cfg.pcc.keepcsd = 'yes';
      cfg.pcc.keepmom = 'yes';
      cfg.pcc.lambda = '5%';
      cfg.frequency = 10;
      cfg.sourcemodel = sourcemodel_eeg;
      cfg.headmodel = headmodel_eeg;
      ft_sourceanalysis(cfg, EEG_freq);
      
      % do dipolefit
      cfg = [];
      cfg.numdipoles = 1;
      % cfg.nonlinear = 'no';
      cfg.sourcemodel = sourcemodel_eeg;
      cfg.headmodel = headmodel_eeg;
      cfg.latency = 0.025;
      ft_dipolefitting(cfg, EEG_tlck_df);
    end
    
  end % function perform_source_analysis


end % main function
