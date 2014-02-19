function inspect_bug1937

% TEST inspect_bug1937
% TEST ft_connectivitysimulation ft_freqanalysis_mvar qsubfeval qsubget

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

%% step 1a

cfg1 = cell(1, nSubjects);

for iSubj = 1:nSubjects
  
  cfg1{iSubj} = [];
  cfg1{iSubj}.order = 5;
  cfg1{iSubj}.toolbox = 'bsmart';
  
  cfg1{iSubj}.inputfile = 'data.mat';
  cfg1{iSubj}.outputfile = sprintf('S%02d_mdata.mat', iSubj);
  
end

% do a sanity check without distributed computing
% mdata = ft_mvaranalysis(cfg{1});

%% step 1b - submit the job to torque

joblist1 = cellfun(@qsubfeval, ...
  repmat({'ft_mvaranalysis'}, [1 nSubjects]), ...
  cfg1, ...
  repmat({'memreq'}, [1 nSubjects]), num2cell(repmat(10*(1024^2), [1 nSubjects])), ...
  repmat({'timreq'}, [1 nSubjects]), num2cell(repmat(30,          [1 nSubjects])), ...
  'UniformOutput', false);

%% step 2a
cfg2 = cell(1, nSubjects);

for iSubj = 1:nSubjects
  
  cfg2{iSubj} = [];
  cfg2{iSubj}.method = 'mvar';
  
  
  cfg2{iSubj}.inputfile = sprintf('S%02d_mdata.mat', iSubj);
  cfg2{iSubj}.outputfile = sprintf('S%02d_mfreq.mat', iSubj);
  
end

% do a sanity check without distributed computing
% mfreq = ft_freqanalysis_mvar(cfg{1});

%% step 2b - submit the job to torque

joblist2 = cellfun(@qsubfeval, ...
  repmat({'ft_freqanalysis_mvar'}, [1 nSubjects]), ...
  cfg2, ...
  repmat({'memreq'}, [1 nSubjects]), num2cell(repmat(10*(1024^2), [1 nSubjects])), ...
  repmat({'timreq'}, [1 nSubjects]), num2cell(repmat(30,          [1 nSubjects])), ...
  repmat({'waitfor'}, [1 nSubjects]), joblist1, ...
  'UniformOutput', false);


%% cleanup

[dum1, options1] = cellfun(@qsubget, joblist1, repmat({'output'}, [1 nSubjects]), repmat({'cell'}, [1 nSubjects]), 'UniformOutput', false);
[dum2, options2] = cellfun(@qsubget, joblist2, repmat({'output'}, [1 nSubjects]), repmat({'cell'}, [1 nSubjects]), 'UniformOutput', false);

cd(prevdir);
rmdir(tmpdir, 's');


