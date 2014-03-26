function inspect_bug1937

% TEST inspect_bug1937
% TEST ft_connectivitysimulation ft_freqanalysis_mvar qsubfeval qsubget

% alebac, 19 Mar 2014: added support "waitfor" for multiple jobs

prevdir = pwd;
tmpdir = tempname(getenv('HOME'));
mkdir(tmpdir);
cd(tmpdir);

nSubjects = 4;


%% generate some data, similar to connectivitextended tutorial
cfg1 = [];
cfg1.ntrials = 10;
cfg1.triallength = 1;
cfg1.fsample = 200;
cfg1.nsignal = 3;
cfg1.method = 'ar';
cfg1.params(:, :, 1) = [ 0.8 0 0;0 0.9 0.5;0.4 0 0.5];
cfg1.params(:, :, 2) = [-0.5 0 0;0 -0.8 0;0 0 -0.2];
cfg1.noisecov = [0.3 0 0;0 1 0;0 0 0.2];

data = ft_connectivitysimulation(cfg1);

% save the data in a temporary directory in your home
save('data.mat', 'data', '-v7.3');

clear data cfg


%% step 1a -  simple preproc

cfg1 = cell(1, nSubjects);

for iSubj = 1:nSubjects
  
  cfg1{iSubj}           = [];
  cfg1{iSubj}.bpfilter  = 'yes';
  cfg1{iSubj}.bpfreq    = [3 30];
  
  cfg1{iSubj}.inputfile = 'data.mat';
  cfg1{iSubj}.outputfile = sprintf('S%02d_data.mat', iSubj);
  
end

% do a sanity check without distributed computing
% data = ft_preprocessing(cfg1{1});


%% step 1b - submit the job to torque

joblist1 = cellfun(@qsubfeval, ...
  repmat({'ft_preprocessing'}, [1 nSubjects]), ...
  cfg1, ...
  repmat({'memreq'}, [1 nSubjects]), num2cell(repmat(30*(1024^2), [1 nSubjects])), ...
  repmat({'timreq'}, [1 nSubjects]), num2cell(repmat(30,          [1 nSubjects])), ...
  'UniformOutput', false);


%% step 2a
cfg2 = cell(1, nSubjects);

for iSubj = 1:nSubjects
  
  cfg2{iSubj}           = [];
  cfg2{iSubj}.method    = 'mtmfft';
  cfg2{iSubj}.taper     = 'hanning';
  
  
  cfg2{iSubj}.inputfile = sprintf('S%02d_data.mat', iSubj);
  cfg2{iSubj}.outputfile = sprintf('S%02d_freq.mat', iSubj);
  
end

% do a sanity check without distributed computing
% freq = ft_freqanalysis(cfg2{1});


%% step 2b - submit the job to torque

joblist2 = cellfun(@qsubfeval, ...
  repmat({'ft_freqanalysis'}, [1 nSubjects]), ...
  cfg2, ...
  repmat({'memreq'},  [1 nSubjects]), num2cell(repmat(10*(1024^2), [1 nSubjects])), ...
  repmat({'timreq'},  [1 nSubjects]), num2cell(repmat(30,          [1 nSubjects])), ...
  repmat({'waitfor'}, [1 nSubjects]), joblist1, ...
  'UniformOutput', false);


%% step 3a - second level analysis

cfg3 = [];
for iSubj = 1:nSubjects
  cfg3.inputfile{iSubj} = sprintf('S%02d_freq.mat', iSubj);
end
cfg3.outputfile = 'grandavg.mat';

% do a sanity check without distributed computing
% grandavg = ft_freqgrandaverage(cfg3)


%% step 3b - submit the job to torque
% should wait for all single subject jobs of step 2 are completed

joblist3 = qsubfeval(@ft_freqgrandaverage, cfg3,'memreq',1*(1024^3),'timreq',60,'waitfor',joblist2);


%% cleanup

[dum1, options1] = cellfun(@qsubget, joblist1, repmat({'output'}, [1 nSubjects]), repmat({'cell'}, [1 nSubjects]), 'UniformOutput', false);
[dum2, options2] = cellfun(@qsubget, joblist2, repmat({'output'}, [1 nSubjects]), repmat({'cell'}, [1 nSubjects]), 'UniformOutput', false);
[dum3, options3] = qsubget(joblist3,'output','cell');

cd(prevdir);
rmdir(tmpdir, 's');


