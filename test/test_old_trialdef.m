function test_old_trialdef

% MEM 1gb
% WALLTIME 00:10:00


% This script tests the implementation of the new representation of trial
% specific information in the data structure.
% The idea is as follows:
% -instead of relying on data.cfg.previous....trl, the trial information
% (when consistent with data representation, so only on the level of raw
% data, and in original sampling rate) should be directly on the main level
% of the structure.
% -this means that raw data should get a field sampleinfo containing the
% first 2 columns of the original trl-matrix.
% -the 3rd column of the trl matrix is the offset column and is represented
% in each trial's individual time axis.
% -user-defined columns > 3 should go to the field trialinfo. this field
% can percolate deeper into the pipeline. yet, once data is not raw
% datatype anymore (or resampled) the sampleinfo field should be removed,
% because consistency with the data is lost.
% -the sampleinfo and trialinfo fields will be updated by subselecting rows
% -sampleinfo should be changed by ft_redefinetrial
% -sampleinfo should be removed by ft_resampledata and any other function
% working on the raw datatype and outputting any other datatype
% -sampleinfo and trialinfo should be concatenated in ft_appenddata
% -sampleinfo and trialinfo should be adjusted by ft_rejectartifact

% needed for the dccnpath function, since we will change directory later on
addpath(fileparts(mfilename('fullpath')));

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'));
headerfile = 'A0132_Aud-Obj-Recognition_20051115_02.res4';
datafile   = 'A0132_Aud-Obj-Recognition_20051115_02.meg4';
%hdr        = ft_read_header(headerfile);

cfg          = [];
cfg.datafile = datafile;
cfg.trl      = [[1001:1000:10001]' [2000:1000:11000]' round(randn(10,1)*100)];
cfg.trl(:,4) = [ones(5,1); ones(5,1)*2];
cfg.continuous = 'yes';
data1          = ft_preprocessing(cfg);
cfg.trl(:,1:2) = cfg.trl(:,1:2) + 11000;
cfg.trl(:,3)   = round(randn(10,1)*100);
data2          = ft_preprocessing(cfg);

% test subselection of trialsclear 
cfg        = [];
cfg.trials = (1:5);
datax      = ft_preprocessing(cfg, data1);

% test some other preprocessing (so nargin==2 for ft_preprocessing)
cfg        = [];
cfg.blc    = 'yes';
datax2     = ft_preprocessing(cfg, data1);

% test checkdata
datay = data1;
datay = rmfield(datay,   'sampleinfo');
datay = rmfield(datay,   'trialinfo');
datay = checkdata(datay, 'hastrialdef', 'yes');

datay.trialinfo(:) = 3;
datay = checkdata(datay, 'hastrialdef', 'yes'); % should give all(datay.trialinfo==3)

datay.trialinfo = [];
datay = checkdata(datay, 'hastrialdef', 'yes'); % should give all(datay.trialinfo==trl(:,4))

datay = data1;
datay = rmfield(datay,   'sampleinfo');
datay = rmfield(datay,   'trialinfo');
datay.cfg = rmfield(datay.cfg, 'trl');
datay = checkdata(datay, 'hastrialdef', 'yes');

% test appenddata
data1b = rmfield(data1, 'sampleinfo');
data1b = rmfield(data1b, 'trialinfo');
%data1b.cfg = rmfield(data1b.cfg, 'trl');
data2b = rmfield(data2, 'sampleinfo');
data2b = rmfield(data2b, 'trialinfo');
%data2b.cfg = rmfield(data2b.cfg, 'trl');
data3 = ft_appenddata([], data1b, data2b);

% test ft_redefinetrial with and without sampleinfo and trialinfo present in data
cfg = [];
cfg.toilim = [-0.5 0.5];
data1bx    = ft_redefinetrial(cfg, data1b);
data1x     = ft_redefinetrial(cfg, data1 );

cfg = [];
cfg.begsample =    1 + round(rand(10,1)*200);
cfg.endsample = 1000 - round(rand(10,1)*200);
data1by    = ft_redefinetrial(cfg, data1b);
data1y     = ft_redefinetrial(cfg, data1 );

cfg = [];
cfg.trl      = data1.cfg.trl;
cfg.trl(:,1) = cfg.trl(:,1)+round(rand(10,1)*20);
cfg.trl(:,2) = cfg.trl(:,2)-round(rand(10,1)*20);
data1bz      = ft_redefinetrial(cfg, data1b);
data1z       = ft_redefinetrial(cfg, data1 );

% test ft_resampledata
cfg = [];
cfg.resamplefs = 300;
cfg.detrend    = 'no';
cfg.blc        = 'yes';
data1rs        = ft_resampledata(cfg, data1);
data1brs       = ft_resampledata(cfg, data1b);
cfg.trials     = [1:5];
data1rs2       = ft_resampledata(cfg, data1b);

% test ft_databrowser
cfg            = [];
cfg.viewmethod = 'butterfly';
cfg.preproc.blc = 'yes';
cfg.channel    = 'MEG';
cfg.continuous = 'yes';
ft_databrowser(cfg, data1b);

% test ft_rejectartifact
cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ctf275/A0132_Aud-Obj-Recognition_20051115_02.ds'));
headerfile = 'A0132_Aud-Obj-Recognition_20051115_02.res4';
datafile   = 'A0132_Aud-Obj-Recognition_20051115_02.meg4';
hdr        = ft_read_header(headerfile);

cfg          = [];
cfg.datafile = datafile;
cfg.trl      = [[1001:1000:10001]' [2000:1000:11000]' round(randn(10,1)*100)];
cfg.trl(:,4) = [ones(5,1); ones(5,1)*2];
cfg.continuous = 'yes';
cfg.artfctdef.type  = 'eog';
cfg.artfctdef.eog.channel = 'MLO11';
cfg.artfctdef.eog.feedback = 'no'; % yes does not work in the test script, see http://bugzilla.fcdonders.nl/show_bug.cgi?id=2840
cfg          = ft_rejectartifact(cfg);

datay = data1;
datay = rmfield(datay,   'sampleinfo');
datay = rmfield(datay,   'trialinfo');
datay2       = ft_rejectartifact(cfg, datay);

% test all other functions using raw data as input (to remove sampleinfo)

% test all high level functions to correctly pass on the trialinfo field
% ft_timelockanalysis
cfg              = [];
cfg.vartrllength = 2;
timelock         = ft_timelockanalysis(cfg, data1);

%ft_denoise_synthetic
cfg          = [];
cfg.gradient = 'G1BR';
dataG1BR     = ft_denoise_synthetic(cfg, data1);
